#htoaa analysis main code

import os
import sys
import json
import glob
from collections import OrderedDict as OD
import time
import tracemalloc
import math
import numpy as np
from copy import deepcopy
#import uproot
import uproot3 as uproot
from auxiliary import met_p4
from auxiliary import top_pT_reweighting
from metfilter import apply_metfilters
np.set_printoptions(threshold=sys.maxsize)
'''
H->aa->4b boosted analysis macro

References:
  * Coffea framework used for TTGamma analysis: https://github.com/nsmith-/TTGamma_LongExercise/blob/FullAnalysis/ttgamma/processor.py
* Coffea installation: /home/siddhesh/anaconda3/envs/ana_htoaa/lib/python3.10/site-packages/coffea
'''

pt = OD([
    ('HT-0To70', [0, 70]),
    ('HT-70To100', [70,100]),
    ('HT-100To200', [100,200]),
    ('HT-200To400', [200,400]),
    ('HT-400To600', [400,600]),
    ('HT-600To800', [600,800]),
    ('HT-800To1200', [800, 1200]),
    ('HT-1200To2500', [1200,2500]),
    ('HT-2500ToInf', [2500]),
])

#import coffea.processor as processor
from coffea import processor, util
from coffea.nanoevents import schemas
from coffea.nanoevents.methods import nanoaod, vector
from coffea.analysis_tools import PackedSelection, Weights
from coffea import hist
import awkward as ak
import uproot


from htoaa_Settings import *#bTagWPs
from htoaa_CommonTools import GetDictFromJsonFile, calculate_lumiScale, setXRootDRedirector, xrdcpFile
from objectselection import ObjectSelection
from genparticle import genparticle

stitchingfile = 'stitchinginfo.json'
with open(stitchingfile) as fSamplesInfo:
    stitchinginfo = json.load(fSamplesInfo)

nEventToReadInBatch = 0.5*10**6 # 2500000 #  1000 # 2500000
nEventsToAnalyze =  -1 #-1 # 1000 # 100000 # -1
sWeighted = "Wtd: "


class HToAATo4bProcessor(processor.ProcessorABC):
    def __init__(self, datasetInfo={}):

        ak.behavior.update(nanoaod.behavior)
        self.datasetInfo = datasetInfo
        #self.isMC = isMC
        self.objectSelector = ObjectSelection(era=self.datasetInfo["era"])
        dataset_axis    = hist.Cat("dataset", "Dataset")
        systematic_axis = hist.Cat("systematic", "Systematic Uncertatinty")

        cutFlow_axis  = hist.Bin("CutFlow",   r"Cuts",            21, -0.5, 20.5)
        genweight_axis  = hist.Bin("genWeight",   r"genWeight",            201, -100.5, 100.5)
        pt_axis       = hist.Bin("Pt",        r"$p_{T}$ [GeV]",   200, 0, 1000)
        mass_axis     = hist.Bin("Mass",      r"$m$ [GeV]",       300, 0, 300)
        mlScore_axis  = hist.Bin("MLScore",   r"ML score",        200, 0., 1.)
        LHE_HT_axis   = hist.Bin("LHE_HT",   r"LHE_HT",        2500, 0., 3000.)
        sXaxis      = 'xAxis'
        sXaxisLabel = 'xAxisLabel'
        sYaxis      = 'yAxis'
        sYaxisLabel = 'yAxisLabel'
        histos = OD([
            ('hCutFlow',                                  {sXaxis: cutFlow_axis,    sXaxisLabel: 'Cuts'}),
            ('hFatJetlen',                                {sXaxis: cutFlow_axis,    sXaxisLabel: 'len(FatJet)'}),
            ('hLeadingFatJetPt',                          {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]"}),
            ('hLeadingFatJetMSoftDrop',                   {sXaxis: mass_axis,       sXaxisLabel: r"m_{soft drop} (leading FatJet) [GeV]"}),
            ('hLeadingFatJetParticleNetMD_Xbb',           {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_Xbb"}),
            ('hLeadingFatJetMass',                        {sXaxis: mass_axis,       sXaxisLabel: r"leading FatJet mass [GeV]"}),
            ('hLeadingFatJetDeepTagMD_bbvsLight',         {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_bbvsLight"}),
            ('LHE_HT_gen',                                {sXaxis: LHE_HT_axis,    sXaxisLabel: r"LHE_HT"}),
            ('genweight',                                  {sXaxis: genweight_axis,    sXaxisLabel: 'genWeight'}),
        ])
        self._accumulator = processor.dict_accumulator({
            'cutflow': processor.defaultdict_accumulator(int)
        })

        for histName, histAttributes in histos.items():
            hXaxis = deepcopy(histAttributes[sXaxis])
            hXaxis.label = histAttributes[sXaxisLabel]
            if sYaxis not in histAttributes.keys():
                # TH1
                self._accumulator.add({
                    histName: hist.Hist(
                        histName,
                        dataset_axis,
                        hXaxis,
                        systematic_axis,
                    )
                })
            else:
                # TH2
                hYaxis = deepcopy(histAttributes[sYaxis])
                hYaxis.label = histAttributes[sYaxisLabel]
                
                self._accumulator.add({
                    histName: hist.Hist(
                        "Counts",
                        dataset_axis,
                        hXaxis,
                        hYaxis,
                        systematic_axis,
                    )
                })
    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        dataset = events.metadata["dataset"] # dataset label
        self.datasetInfo[dataset]['isSignal'] = False
        self.datasetInfo[dataset]['isQCD'] = False

        events = apply_metfilters(self.datasetInfo[dataset]['isMC'], events)
        if self.datasetInfo[dataset]['isMC']:
            self.datasetInfo[dataset]['isSignal'] = True if "HToAATo4B" in dataset else False
            self.datasetInfo[dataset]['isQCD']    = True if "QCD"       in dataset else False
            
            output = self.accumulator.identity()
            systematics_shift = [None]
            for _syst in systematics_shift:
                output += self.process_shift(events, _syst)
        else:
            output = self.process_shift(events, None)

        return output
    
    def process_shift(self, events, shift_syst=None):

        output = self.accumulator.identity()
        dataset = events.metadata["dataset"] # dataset label
        ##################
        # OBJECT SELECTION
        ##################


        # Gen-level selection ---------------------------------------------------------------------
        genHiggs = None
        genHT    = None
        if self.datasetInfo[dataset]['isMC'] and self.datasetInfo[dataset]['isSignal']: 
            genHiggs  = self.objectSelector.selectGenHiggs(events)        
            genHT     = self.objectSelector.GenHT(events)
            # m(bbar from A) and m(4b from HToAA) ----------

            genACollection = self.objectSelector.selectGenABoson(events)
            genA_First  = genACollection[:, 0]
            genA_Second = genACollection[:, 1]
            
            idxGenA_sortByMass = ak.argsort(genACollection.mass, axis=-1, ascending=False)
            
            genBBar_pairs_all = ak.argcombinations(events.GenPart, 2, fields=['b', 'bbar'])
            genBBar_pairs = genBBar_pairs_all[(
                (abs(events.GenPart[genBBar_pairs_all['b'   ]].pdgId) == 5) &
                (abs(events.GenPart[genBBar_pairs_all['bbar']].pdgId) == 5) &
                ((events.GenPart[genBBar_pairs_all['b']].pdgId) == (-1*events.GenPart[genBBar_pairs_all['bbar']].pdgId)  ) &
                (events.GenPart[genBBar_pairs_all['b']].genPartIdxMother == events.GenPart[genBBar_pairs_all['bbar']].genPartIdxMother) &
                (events.GenPart[ events.GenPart[genBBar_pairs_all['b'   ]].genPartIdxMother ].pdgId == 36) &
                (events.GenPart[ events.GenPart[genBBar_pairs_all['bbar']].genPartIdxMother ].pdgId == 36) &
                (events.GenPart[genBBar_pairs_all['bbar']].genPartIdxMother != -1) &
                (events.GenPart[genBBar_pairs_all['b']].genPartIdxMother != -1)
            )]
            # LorentVector of GenB quarks from HToAATo4b
            nEvents_11 = ak.num(events.GenPart[genBBar_pairs['b']][:, 0].pt, axis=0)

            # https://coffeateam.github.io/coffea/modules/coffea.nanoevents.methods.vector.html
            LVGenB_0 = genparticle(events, genBBar_pairs['b'], 0)
            LVGenBbar_0 = genparticle(events, genBBar_pairs['bbar'], 0)
            LVGenB_1 = genparticle(events, genBBar_pairs['b'], 1)
            LVGenBbar_1 = genparticle(events, genBBar_pairs['bbar'], 1) 
        # QCD MC ----------------------------------------------
        if self.datasetInfo[dataset]['isMC'] and self.datasetInfo[dataset]['isQCD'] :
            genBQuarks_QCD = events.GenPart[(
                (abs(events.GenPart.pdgId) == 5)
            )]

        # Reco-level -----------------------------------------------------------------------------------
        
        ##################
        # EVENT VARIABLES
        ##################

        HLT_AK8PFJet330_name = "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4" 
        
        # sel_names_all = dict of {"selection name" : [list of different cuts]}; for cut-flow table
        sel_names_all = OD([
            ("SR",                    [
                "nPV",
                "FatJet",
                "lep",
                "met"
                #"ak4Jet"
            ]),
        ])
        # reconstruction level cuts for cut-flow table. Order of cuts is IMPORTANT
        cuts_reco = ["dR_LeadingFatJet_GenB_0p8"] + sel_names_all["SR"] #.copy()
        selection = PackedSelection()
        FatJet = self.objectSelector.selectFatJets(events)
        ak4Jet = self.objectSelector.selectak4Jets(events)
        ele = self.objectSelector.selectElectrons(events)
        mu = self.objectSelector.selectMuons(events)
        lep = ak.with_name(
                ak.concatenate([mu, ele], axis=1), "PtEtaPhiMCandidate"
            )
        lep_et = np.sqrt(lep.mass * lep.mass + lep.px * lep.px + lep.py * lep.py)
        m = met_p4(events)
        #nu_et = np.sqrt(m.px * m.px + m.py*m.py)
        mt_W = lep.mass * lep.mass +  2* (lep.pt*m.pt -(lep.px *m.px + lep.py * m.py))
        # create a PackedSelection object
        # this will help us later in composing the boolean selections easily
        if self.datasetInfo[dataset]['isMC'] and self.datasetInfo[dataset]['isSignal']:
            sel_names_GEN = ["1GenHiggs", "2GenA", "2GenAToBBbarPairs", "dR_GenH_GenB_0p8"]
            FatJet = FatJet[(FatJet.delta_r(LVGenB_0) <0.8)
                            & (FatJet.delta_r(LVGenB_1) <0.8)
                            & (FatJet.delta_r(LVGenBbar_0) < 0.8)
                            & (FatJet.delta_r(LVGenBbar_1) < 0.8)
            ]
            selection.add("1GenHiggs", ak.num(genHiggs) == 1)
            selection.add("2GenA", ak.num(genACollection) == 2)
            selection.add("2GenAToBBbarPairs", ak.num(genBBar_pairs) == 2)

        FatJet = FatJet[ak.all(FatJet.metric_table(lep)>0.8, axis=-1)]
        if "nPV" in sel_names_all["SR"]:
            selection.add("nPV", events.PV.npvsGood >= 1)
        if "FatJet" in sel_names_all["SR"]:
            selection.add(
                "FatJet",
                ak.num(FatJet, axis=1) >=1
            )
        if 'lep' in sel_names_all["SR"]:
            selection.add("lep", ak.num(lep) == 1)
        if 'met' in sel_names_all["SR"]:
            selection.add("met", events.MET.pt > 20)
        if 'ak4Jet' in sel_names_all["SR"]:
            selection.add("ak4Jet", ak.num(ak4Jet) < 1)
        # useful debugger for selection efficiency
        sel_SR           = selection.all(* sel_names_all["SR"])
        ################
        # EVENT WEIGHTS
        ################
        # create a processor Weights object, with the same length as the number of events in the chunk
        weights     = Weights(len(events))
        weights_gen = Weights(len(events))
        genweight_unique = set(events.genWeight)
        count = {}
        for gu in genweight_unique:
            count[gu] = ak.sum(ak.where(events.genWeight==gu, 1,0))
        count = sorted([(k,v) for k,v in count.items()], key=lambda kv: kv[1], reverse=True)
        events.genWeight = ak.where(events.genWeight >3*count[0][0], 3*count[0][0], events.genWeight)
        stitch = np.full(len(events), self.datasetInfo[dataset]["lumiScale"])
        for (k,v) in pt.items():
            if len(v) > 1:
                stitch = ak.where((events.LHE.HT>=v[0])&(events.LHE.HT<v[1]), 59.83*1000*stitchinginfo[k]['xs']/stitchinginfo[k]['nevent'], stitch)
            else:
                stitch = ak.where((events.LHE.HT>=v[0]), 59.83*1000*stitchinginfo[k]['xs']/stitchinginfo[k]['nevent'], stitch)
        if self.datasetInfo[dataset]["isMC"]:
            if not self.datasetInfo[dataset]["stitching"]:
                weights.add(
                    "lumiWeight",
                    weight=np.full(len(events), self.datasetInfo[dataset]["lumiScale"])
                )
                weights_gen.add(
                    "lumiWeight",
                    weight=np.full(len(events), self.datasetInfo[dataset]["lumiScale"])
                )
            else:
                weights.add(
                    "lumiWeight",
                    weight=stitch#np.full(len(events), self.datasetInfo[dataset]["lumiScale"])
                )
                weights_gen.add(
                    "lumiWeight",
                    weight=stitch#np.full(len(events), self.datasetInfo[dataset]["lumiScale"])
                )
            weights.add(
                "genWeight",
                weight=events.genWeight#np.copysign(np.ones(len(events)), events.genWeight)
            )
            weights_gen.add(
                "genWeight",
                weight=events.genWeight#np.copysign(np.ones(len(events)), events.genWeight)
            )
            weights.add(
                "btagWeight",
                weight=(events.btagWeight.DeepCSVB)
            )
            if dataset.startswith('TT'):
                topt_reweiging = ak.to_numpy(top_pT_reweighting(events.GenPart)).squeeze()
                weights.add(
                    "topt_reweiging",
                    weight=topt_reweiging,
                    weightUp=topt_reweiging**2,  # w**2
                    weightDown=ak.ones_like(topt_reweiging)  # w = 1
                )
            
        ###################
        # FILL HISTOGRAMS
        ###################
        systList = []
        if self.datasetInfo[dataset]['isMC']:
            if shift_syst is None:
                systList = [
                    "central"
                ]
                if dataset.startswith('TT'):
                   systList.extend(["topt_reweigingUp", "topt_reweigingDown"])
            else:
                systList = [shift_syst]
        else:
            systList = ["noweight"]
        #systList = ["noweight"]
            
        output['cutflow']['all events'] += len(events)
        output['cutflow'][sWeighted+'all events'] += weights.weight().sum()
        for iSelection in sel_names_all.keys():
            iName = f"{iSelection}: {sel_names_all[iSelection]}"
            sel_i = selection.all(* sel_names_all[iSelection])
            output['cutflow'][iName] += sel_i.sum()
            output['cutflow'][sWeighted+iName] +=  weights.weight()[sel_i].sum()
        allcuts = sel_names_all["SR"]
        for syst in systList:

            # find the event weight to be used when filling the histograms
            weightSyst = syst
            weightSyst_gen = None
            # in the case of 'central', or the jet energy systematics, no weight systematic variation is used (weightSyst=None)
            if syst in ["central", "JERUp", "JERDown", "JESUp", "JESDown"]:
                weightSyst = None

            ones_list = np.ones(len(events))
            
            if syst == "noweight":
                evtWeight = np.ones(len(events))
            else:
                evtWeight     = weights.weight(weightSyst)
                evtWeight_gen = weights_gen.weight(weightSyst_gen)
            # all events
            cuts = []
            for ibin, cut in enumerate(allcuts):
                cuts.append(cut)
                n = selection.all(*cuts).sum()
                n = np.ones(n)
                output['hCutFlow'].fill(
                    dataset=dataset,
                    CutFlow=(n * ibin),
                    systematic=syst
                )
            output['hLeadingFatJetPt'].fill(
                dataset=dataset,
                Pt=(ak.flatten(FatJet.pt[sel_SR][:, 0:1])),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['LHE_HT_gen'].fill(
                dataset=dataset,
                LHE_HT=(events.LHE.HT),
                systematic=syst,
                weight=evtWeight_gen
            )
            output['genweight'].fill(
                dataset=dataset,
                genWeight=(events.genWeight),
                systematic=syst
            )
            output['hLeadingFatJetMSoftDrop'].fill(
                dataset=dataset,
                Mass=(ak.flatten(FatJet.msoftdrop[sel_SR][:, 0:1])),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hFatJetlen'].fill(
                dataset=dataset,
                CutFlow=(ak.num(FatJet[sel_SR])),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            fatjet_sortedpnmd = FatJet[ak.argsort(FatJet.particleNetMD_Xbb, ascending=False)]
            output['hLeadingFatJetParticleNetMD_Xbb'].fill(
                dataset=dataset,
                MLScore=(ak.flatten(fatjet_sortedpnmd.particleNetMD_Xbb[sel_SR][:, 0:1])), 
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            fatjet_sortedbbvs = FatJet[ak.argsort(FatJet.deepTagMD_bbvsLight, ascending=False)]
            output['hLeadingFatJetDeepTagMD_bbvsLight'].fill(
                dataset=dataset,
                MLScore=(ak.flatten(fatjet_sortedbbvs.deepTagMD_bbvsLight[sel_SR][:, 0:1])),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetMass'].fill(
                dataset=dataset,
                Mass=(ak.flatten(fatjet_sortedbbvs.mass[sel_SR][:, 0:1])),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
        return output


    def postprocess(self, accumulator):
        #pass
        return accumulator

if __name__ == '__main__':
    print("htoaa_Analysis:: main: {}".format(sys.argv)); sys.stdout.flush()

    if len(sys.argv) != 2:
        print("htoaa_Analysis:: Command-line config file missing.. \t **** ERROR **** \n")

    sConfig = sys.argv[1]

    config = GetDictFromJsonFile(sConfig)
    print("Config {}: \n{}".format(sConfig, json.dumps(config, indent=4)))

    lumiScale = 1
    sInputFiles         = config["inputFiles"]
    sOutputFile         = config["outputFile"]
    sample_category     = config['sampleCategory']
    process_name        = config['process_name']
    isMC                = config["isMC"]
    era                 = config['era']
    if isMC:
        luminosity          = Luminosities[era][0]
        sample_crossSection = config["crossSection"]
        sample_nEvents      = config["nEvents"]
        sample_sumEvents    = config["sumEvents"] if config["sumEvents"] != -1 else sample_nEvents
        if sample_sumEvents == -1: sample_sumEvents = 1 # Case when sumEvents is not calculated
        lumiScale = calculate_lumiScale(luminosity=luminosity, crossSection=sample_crossSection, sumEvents=sample_sumEvents)

    sInputFiles_toUse = []
    for sInputFile in sInputFiles:
        if "*" in sInputFile:
            sInputFiles_toUse.extend( glob.glob( sInputFile ) )
        elif 'eos' not in sInputFile:
            sInputFile = setXRootDRedirector(sInputFile)
            sFileLocal = f'/tmp/snandan/inputFiles/{process_name}/{os.path.basename(sInputFile)}'
            if xrdcpFile(sInputFile, sFileLocal, nTry = 3):
                sInputFiles_toUse.append(sFileLocal)
            else:
                print(f"Ip file {sInputFile} failed to download \t **** ERROR ****")
                exit(1)
        else:
            sInputFiles_toUse.append( sInputFile )
    sInputFiles = sInputFiles_toUse
    sys.stdout.flush()
    startTime = time.time()
    tracemalloc.start()
    chunksize = nEventToReadInBatch
    maxchunks = None if nEventsToAnalyze == -1 else int(nEventsToAnalyze/nEventToReadInBatch)
    print(f"nEventsToAnalyze: {nEventsToAnalyze},  nEventToReadInBatch: {nEventToReadInBatch}, chunksize: {chunksize},  maxchunks: {maxchunks}")
    run = processor.Runner(
        #executor=executor,
        executor=processor.FuturesExecutor(workers=4),
        schema=schemas.NanoAODSchema,
        savemetrics=True,
        chunksize=chunksize,  #3 ** 20,  ## Governs the number of times LeptonJetProcessor "process" is called
        maxchunks=maxchunks
    )

    output, metrics = run(
        fileset={sample_category: sInputFiles},
        treename="Events",
        processor_instance=HToAATo4bProcessor(
            datasetInfo={
                "era": era,
                sample_category: {"isMC": isMC, "lumiScale": lumiScale, "stitching":True}
            }
        )
    )

    print(f"metrics: {metrics}")
    
    if 'cutflow' in output.keys():
        print("Cutflow::")
        for key in output['cutflow'].keys():
            if key.startswith(sWeighted): continue
            print("%10f\t%10d\t%s" % (output['cutflow'][sWeighted+key], output['cutflow'][key], key))
    
    if sOutputFile is not None:
        if not sOutputFile.endswith('.root'): sOutputFile += '.root'
        sDir1 = 'evt/%s' % (sample_category)
        with uproot.recreate(sOutputFile) as fOut:
            for key, value in output.items():
                if not isinstance(value, hist.Hist): continue
                for _dataset in value.axis('dataset').identifiers():
                    for _syst in value.axis('systematic').identifiers():
                        h1 = value.integrate('dataset',_dataset).integrate('systematic',_syst).to_hist()
                        if 'central' not in _syst.name:
                            histname = f'{key}_{_syst}'
                        else:
                            histname = key
                        if 'gen' not in value.label:
                            fOut[f'{sDir1}/{histname}'] = h1
                        else:
                            if 'central' not in _syst.name: continue
                            fOut[f'{sDir1}/gen/{histname}'] = h1
    current_memory, peak_memory = tracemalloc.get_traced_memory() # https://medium.com/survata-engineering-blog/monitoring-memory-usage-of-a-running-python-program-49f027e3d1ba
    print(f"\n\nMemory usage:: current {current_memory / 10**6}MB;  peak {peak_memory / 10**6}MB")

    endTime = time.time()
    totalTime = endTime - startTime
    totalTime_hr  = int(totalTime/60/60)
    totalTime_min = totalTime - float(totalTime_hr * 60)
    totalTime_min = int(totalTime_min/60)
    totalTime_sec = totalTime - float(totalTime_hr * 60*60) - float(totalTime_min * 60)
    os.system(f'rm -r /tmp/snandan/inputFiles/{process_name}')
    print(f"Total run time: {totalTime_hr}h {totalTime_min}m {totalTime_sec}s = {totalTime}sec ")
