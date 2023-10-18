#https://cms.cern.ch/iCMS/analysisadmin/cadilines?line=HIG-18-026&tp=an&id=2158&ancode=HIG-18-026
#htoaa analysis main code
#/afs/cern.ch/work/s/snandan/anaconda3/envs/myEnv/lib/python3.10/site-packages/boost_histogram/_internal/hist.py

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
from auxiliary import selectRunLuminosityBlock
from metfilter import apply_metfilters
np.set_printoptions(threshold=sys.maxsize)
import hist
from histograms import histograms, clip_value
from weights import event_weights#calculate_weight, transfer_factor
'''
H->aa->4b boosted analysis macro

References:
  * Coffea framework used for TTGamma analysis: https://github.com/nsmith-/TTGamma_LongExercise/blob/FullAnalysis/ttgamma/processor.py
* Coffea installation: /home/siddhesh/anaconda3/envs/ana_htoaa/lib/python3.10/site-packages/coffea
'''
#import coffea.processor as processor
from coffea import processor, util
from coffea.nanoevents import schemas
from coffea.nanoevents.methods import nanoaod, vector
from coffea.analysis_tools import PackedSelection, Weights
#from coffea import hist
import awkward as ak
import uproot


from htoaa_Settings import *#bTagWPs
from htoaa_CommonTools import GetDictFromJsonFile, calculate_lumiScale, setXRootDRedirector, xrdcpFile
from objectselection import ObjectSelection
from genparticle import genparticle

nEventToReadInBatch = 0.5*10**6 # 2500000 #  1000 # 2500000
nEventsToAnalyze =  -1 #-1 # 1000 # 100000 # -1
sWeighted = "Wtd: "


class HToAATo4bProcessor(processor.ProcessorABC):
    def __init__(self, datasetInfo={}):

        ak.behavior.update(nanoaod.behavior)
        self.config = datasetInfo
        self.objectSelector = ObjectSelection(era=self.config["era"])

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        dataset = events.metadata["dataset"] # dataset label
        if not self.config['isMC']:
            with open('Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt') as f:
                rl = json.load(f)
            events = events[selectRunLuminosityBlock(rl, events.run, events.luminosityBlock)]
        events = apply_metfilters(self.config['isMC'], events)
        # Gen-level selection ---------------------------------------------------------------------
        genHiggs = None
        genHT    = None
        
        if self.config['isSignal']:
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
        if self.config['isQCD'] :
            genBQuarks_QCD = events.GenPart[(
                (abs(events.GenPart.pdgId) == 5)
            )]

        # Reco-level -----------------------------------------------------------------------------------
        
        ##################
        # EVENT VARIABLES
        ##################

        FatJet = self.objectSelector.selectFatJets(events)
        ak4Jet = self.objectSelector.selectak4Jets(events)
        preselele = self.objectSelector.selectElectrons(events)
        preselmu = self.objectSelector.selectMuons(events)
        elesize = ak.sum(ak.num(preselele))
        musize = ak.sum(ak.num(preselmu))
        is_triggered_1ele = events.HLT.Ele32_WPTight_Gsf
        is_triggered_1mu = events.HLT.IsoMu24 | events.HLT.IsoMu27
        sel_triggered_1ele = (self.config['use_triggers_1e']) & is_triggered_1ele
        sel_triggered_1mu = (self.config['use_triggers_1mu']) & is_triggered_1mu
        preselele = ak.with_field(preselele, (preselele.mvaTTH >=0.3), "tight_ele")
        preselmu = ak.with_field(preselmu, (preselmu.mvaTTH >=0.5), "tight_mu")
        preselele = ak.with_field(preselele, (preselele.mvaTTH <0.3), "fake_ele")
        preselmu = ak.with_field(preselmu, (preselmu.mvaTTH <0.5), "fake_mu")

        #https://github.com/HEP-KBFI/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L94 
        trig_obj_ele = (abs(events.TrigObj.id) == 11) & ( (events.TrigObj.filterBits & 2) == 2 )
        trig_obj_mu = (abs(events.TrigObj.id) == 13) & ( (events.TrigObj.filterBits & 8) == 8 )

        if self.config['isSignal']:
            sel_names_GEN = ["1GenHiggs", "2GenA", "2GenAToBBbarPairs", "dR_GenH_GenB_0p8"]
            FatJet = FatJet[(FatJet.delta_r(LVGenB_0) <0.8)
                            & (FatJet.delta_r(LVGenB_1) <0.8)
                            & (FatJet.delta_r(LVGenBbar_0) < 0.8)
                            & (FatJet.delta_r(LVGenBbar_1) < 0.8)
            ]
            '''selection.add("1GenHiggs", ak.num(genHiggs) == 1)
            selection.add("2GenA", ak.num(genACollection) == 2)
            selection.add("2GenAToBBbarPairs", ak.num(genBBar_pairs) == 2)'''
        FatJet = ak.with_field(FatJet, FatJet.msoftdrop >= 90, "msoftdrop_GE_90")
        FatJet = ak.with_field(FatJet, FatJet.msoftdrop < 90, "msoftdrop_L_90")
        evt_weights = event_weights(events, self.config)
        weights = evt_weights.calculate_weight()

        systList = ["central"]
        '''if self.config['isMC']:
        if dataset.startswith('TT'):
        systList.extend(["topt_reweigingUp", "topt_reweigingDown"])'''

        output = histograms(systList, self.config['msoftdrop'])
        output['cutflow'] = {}
        output['cutflow_weighted'] = {}
        #for lepton in self.config['lepton_selection']:#['tight_lep', 'fake_lepton']:
        if self.config['lepton_selection'] == 'tight_lep':
            selele = preselele[preselele.tight_ele]
            selmu = preselmu[preselmu.tight_mu]
            sellep_all = ak.with_name(
                ak.concatenate([selmu, selele], axis=1), "PtEtaPhiMCandidate"
            )
        else:
            if 1:
                sellep_all = preselele[preselele.fake_ele]
                sellep_tight = preselele[preselele.tight_ele]
            else:
                sellep_all = preselmu[preselmu.fake_mu]
                sellep_tight = preselmu[preselmu.tight_mu]
        sellep_all = sellep_all[ak.argsort(sellep_all.pt, axis=-1, ascending=False)]
        sellep = sellep_all[:,0:1]
        sel_ele = abs(sellep.pdgId) == 11
        sel_mu = abs(sellep.pdgId) == 13
        trigObj_matched_ele = ak.any(events.TrigObj[trig_obj_ele].metric_table(sellep[sel_ele]) < 0.4, axis=-1)
        trigObj_matched_mu = ak.any(events.TrigObj[trig_obj_mu].metric_table(sellep[sel_mu]) < 0.4, axis=-1)
        for msoftdrop in self.config['msoftdrop']:#['msoftdrop_GE_90', 'msoftdrop_L_90']:
            if msoftdrop == 'msoftdrop_GE_90':
                selFatJet = FatJet[(FatJet.msoftdrop > 90)
                & (ak.all(FatJet.metric_table(sellep)>0.8, axis=-1))
                ]
            else:
                assert(msoftdrop == 'msoftdrop_L_90')
                selFatJet = FatJet[(FatJet.msoftdrop < 90)
                    & (ak.all(FatJet.metric_table(sellep)>0.8, axis=-1))
                ]
            ak4jet_cleaned_wrt_lepton = ak4Jet[ak.all(ak4Jet.metric_table(sellep)>0.4, axis=-1)]
            ak4jet_cleaned = ak4jet_cleaned_wrt_lepton[
                ak.all(ak4jet_cleaned_wrt_lepton.metric_table(selFatJet[:,0:1])>0.8, axis=-1)
            ]
            nbtag_medium = ak.num(ak4jet_cleaned[ak4jet_cleaned.btagDeepFlavB > 0.2770])
            ak4jet_notcleaned = ak4Jet[
                ak.all(ak4Jet.metric_table(selFatJet[:,0:1])<=0.8, axis=-1)
            ]
            nbtag_notcleaned_medium = ak.num(ak4jet_notcleaned[ak4jet_notcleaned.btagDeepFlavB > 0.2770])
            ht = ak.sum(ak4jet_cleaned.pt, axis=1)
            
            selection = PackedSelection()
            selection.add("nPV", events.PV.npvsGood >= 1)
            if self.config['lepton_selection'] == 'tight_lep':
                selection.add( self.config['lepton_selection'], (ak.num(sellep_all) == 1) & (ak.num(sellep) == 1))
            else:
                selection.add( self.config['lepton_selection'], (ak.num(sellep_all) > 0)\
                    & (ak.num(sellep_tight) == 0)\
                )
            selection.add('nbtag_medium==0', nbtag_medium == 0)
            if msoftdrop == 'msoftdrop_GE_90':
                selection.add(msoftdrop, ak.num(selFatJet) > 0)
            else:
                selection.add(msoftdrop, (ak.num(selFatJet) > 0)\
                  & (ak.num(FatJet[FatJet.msoftdrop > 90]) == 0))
            selection.add("met", events.MET.pt >= 30)
            selection.add("triggered", sel_triggered_1ele | sel_triggered_1mu)
            selection.add("trigObj_matched_lep",\
               (ak.any(trigObj_matched_ele, axis=-1) & (sel_triggered_1ele))\
                | (ak.any(trigObj_matched_mu, axis=-1) & sel_triggered_1mu ))
            if not self.config['isMC']:
                selection.add("dataset_check", (sel_triggered_1ele & is_triggered_1mu)==0)
            f0 = np.ones(len(events))
            if not self.config['isSignal']\
               and self.config['lepton_selection'] == 'tight_lep'\
               and msoftdrop == 'msoftdrop_L_90'\
               and len(events):
                f0 = evt_weights.transfer_factor(ht, f0)

            allcuts = selection.names#commoncuts + [lepton] + [msoftdrop]
            sel_reg = selection.all(*allcuts)
            ###################
            # FILL HISTOGRAMS
            ###################
            for syst in systList:
                weightSyst = syst
                if syst in ["central", "JERUp", "JERDown", "JESUp", "JESDown"]:
                    weightSyst = None
                if not self.config['isMC']:
                    evtWeight = f0 if self.config['lepton_selection'] == 'tight_lep'\
                                and msoftdrop == 'msoftdrop_L_90'\
                                else np.ones(len(events))
                else:
                    if self.config['lepton_selection'] == 'tight_lep'\
                       and msoftdrop == 'msoftdrop_L_90':
                        evtWeight = weights.weight(weightSyst)*f0
                    else:
                        evtWeight = weights.weight(weightSyst)
                    evtWeight_gen = weights.partial_weight(['lumiWeight', 'genWeight'])
                    
                cuts = []
                if syst == 'central':
                    for ibin, cut in enumerate(allcuts):
                        cuts.append(cut)
                        n = selection.all(*cuts).sum()
                        if msoftdrop not in output['cutflow'].keys():
                            output['cutflow'][msoftdrop] = {}
                            output['cutflow_weighted'][msoftdrop] = {}
                        output['cutflow'][msoftdrop][cut] = n
                        output['cutflow_weighted'][msoftdrop][cut] = np.sum(evtWeight[selection.all(*cuts)])
                        print('cut: ', cut, '\t', n)
                        n = np.ones(n)
                        output['hCutFlow'].fill(
                            jet=msoftdrop,
                            cutflow = (n*ibin),
                            systematic='central',
                            weight=evtWeight[selection.all(*cuts)]
                        )
                output['hLeadingFatJetPt'].fill(
                    jet=msoftdrop,
                    hLeadingFatJetPt=clip_value(
                        ak.flatten(selFatJet.pt[sel_reg][:, 0:1]),
                        min(output['hLeadingFatJetPt'].axes[2].edges),
                        max(output['hLeadingFatJetPt'].axes[2].edges)
                    ),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                output['hLeadingFatJetMSoftDrop'].fill(
                    jet=msoftdrop,
                    hLeadingFatJetMSoftDrop=clip_value(
                        ak.flatten(selFatJet.msoftdrop[sel_reg][:, 0:1]),
                        min(output['hLeadingFatJetMSoftDrop'].axes[2].edges),
                        max(output['hLeadingFatJetMSoftDrop'].axes[2].edges)
                        ),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                output['HT'].fill(
                    jet=msoftdrop,
                    HT=clip_value(
                        ht[sel_reg],
                        min(output['HT'].axes[2].edges),
                        max(output['HT'].axes[2].edges)
                    ),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                if self.config['isMC'] and syst == 'central' and 'LHE' in events.fields:
                    output['LHE_HT_gen'].fill(
                        jet=msoftdrop,
                        LHE_HT_gen=clip_value(
                            events.LHE.HT,
                            min(output['LHE_HT_gen'].axes[2].edges),
                            max(output['LHE_HT_gen'].axes[2].edges)
                        ),
                        systematic=syst,
                        weight=evtWeight_gen
                    )
                    output['LHE_Njet_gen'].fill(
                        jet=msoftdrop,
                        LHE_Njet_gen=(events.LHE.Njets),
                        systematic=syst,
                        weight=evtWeight_gen
                    )
                    output['genweight'].fill(
                        jet=msoftdrop,
                        genWeight=(events.genWeight),
                        systematic=syst
                    )
                fatjet_sortedpnmd = selFatJet[ak.argsort(selFatJet.particleNetMD_Xbb, ascending=False)][:,0:1]
                deltar = fatjet_sortedpnmd.subjets[:,:,0:1].delta_r(fatjet_sortedpnmd.subjets[:,:,1:2])
                output['hLeadingFatJetParticleNetMD_Xbb'].fill(
                    jet=msoftdrop,
                    hLeadingFatJetParticleNetMD_Xbb=(ak.flatten(fatjet_sortedpnmd.particleNetMD_Xbb[sel_reg][:, 0:1])),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                if syst == 'central':
                    cut_womet = ['nPV', msoftdrop]#, 'triggered']
                    ele = events.Electron#[events.Electron.mvaTTH > -0.6]
                    mu = events.Muon#[events.Muon.mvaTTH > -0.6]
                    mva_lep = ak.with_name(ak.concatenate([ele,mu], axis=1), "PtEtaPhiMCandidate")
                    mva_lep = mva_lep[ak.argsort(mva_lep.pt, ascending=False)][:,0:1]
                    output['hmet'].fill(
                        jet=msoftdrop,
                        met=np.minimum(
                            events.MET.pt[selection.all(*cut_womet)],
                            max(output['hmet'].axes[2].edges)-1
                        ),
                        systematic=syst,
                        weight=evtWeight[selection.all(*cut_womet)]
                    )
                    selection_nocut = PackedSelection()
                    selection_nocut.add("lep", ak.num(mva_lep) >= 1)
                    selection_nocut.add("fatjet", ak.num(events.FatJet) >= 1)
                    output['hmet_vs_lepmva'].fill(
                        jet=msoftdrop,
                        met=events.MET.pt[selection_nocut.all(*["lep"])],
                        lepmva=(ak.flatten(mva_lep.mvaTTH[selection_nocut.all(*["lep"])])),
                        systematic=syst,
                        weight=evtWeight[selection_nocut.all(*['lep'])]
                    )
                    output['hmet_vs_lepiso'].fill(
                        jet=msoftdrop,
                        met=events.MET.pt[selection_nocut.all(*["lep"])],
                        lepiso=(ak.flatten(mva_lep.miniPFRelIso_all[selection_nocut.all(*["lep"])])),
                        systematic=syst,
                        weight=evtWeight[selection_nocut.all(*['lep'])]
                    )
                    output['hsoftdrop_vs_lepiso'].fill(
                        jet=msoftdrop,
                        softdrop=(ak.flatten(events.FatJet[:,0:1].msoftdrop[selection_nocut.all(*["lep", 'fatjet'])])),
                        lepiso=(ak.flatten(mva_lep.miniPFRelIso_all[selection_nocut.all(*["lep",'fatjet'])])),
                        systematic=syst,
                        weight=evtWeight[selection_nocut.all(*['lep','fatjet'])]
                    )
                    output['hsoftdrop_vs_lepmva'].fill(
                        jet=msoftdrop,
                        softdrop=(ak.flatten(events.FatJet[:,0:1].msoftdrop[selection_nocut.all(*["lep", 'fatjet'])])),
                        lepmva=(ak.flatten(mva_lep.mvaTTH[selection_nocut.all(*["lep",'fatjet'])])),
                        systematic=syst,
                        weight=evtWeight[selection_nocut.all(*['lep','fatjet'])]
                    )
                    output['leptonsize'].fill(
                        jet=msoftdrop,
                        leptonsize=(np.concatenate((np.ones(elesize), np.ones(musize)*2))),
                        systematic=syst
                    )
                    output['ak4jetsize'].fill(
                        jet=msoftdrop,
                        ak4jetsize=(ak.num(ak4jet_cleaned)[sel_reg]),
                        systematic=syst,
                        weight=evtWeight[sel_reg]
                    )
                    output['ak4jetsize_notcleaned'].fill(
                        jet=msoftdrop,
                        ak4jetsize_notcleaned=(np.minimum(
                            ak.num(ak4jet_notcleaned)[sel_reg],
                            max(output['ak4jetsize_notcleaned'].axes[2].edges)-1
                        )
                    ),
                        systematic=syst,
                        weight=evtWeight[sel_reg]
                    )
                    output['ak4jetsize_mediumbtag'].fill(
                        jet=msoftdrop,
                        ak4jetsize_mediumbtag=(nbtag_medium[sel_reg]),
                        systematic=syst,
                        weight=evtWeight[sel_reg]
                    )
                    output['ak4jetsize_notcleaned_mediumbtag'].fill(
                        jet=msoftdrop,
                        ak4jetsize_notcleaned_mediumbtag=(nbtag_notcleaned_medium[sel_reg]),
                        systematic=syst,
                        weight=evtWeight[sel_reg]
                    )
                    output['tau1'].fill(
                        jet=msoftdrop,
                        tau1=ak.flatten(fatjet_sortedpnmd.tau1[sel_reg][:,0:1]),
                        systematic=syst,
                        weight=evtWeight[sel_reg]
                    )
                    output['tau2'].fill(
                        jet=msoftdrop,
                        tau2=ak.flatten(fatjet_sortedpnmd.tau2[sel_reg][:,0:1]),
                        systematic=syst,
                        weight=evtWeight[sel_reg]
                    )
                    output['tau3'].fill(
                        jet=msoftdrop,
                        tau3=ak.flatten(fatjet_sortedpnmd.tau3[sel_reg][:,0:1]),
                        systematic=syst,
                        weight=evtWeight[sel_reg]
                    )
                    output['tau4_2'].fill(
                        jet=msoftdrop,
                        tau4_2=ak.flatten((fatjet_sortedpnmd.tau4/fatjet_sortedpnmd.tau2)[sel_reg][:,0:1]),
                        systematic=syst,
                        weight=evtWeight[sel_reg]
                    )
                    output['tau4_3'].fill(
                        jet=msoftdrop,
                        tau4_3=ak.flatten((fatjet_sortedpnmd.tau4/fatjet_sortedpnmd.tau3)[sel_reg][:,0:1]),
                        systematic=syst,
                        weight=evtWeight[sel_reg]
                    )
                    output['nConstituents'].fill(
                        jet=msoftdrop,
                        nConstituents=ak.flatten(fatjet_sortedpnmd.nConstituents[sel_reg][:,0:1]),
                        systematic=syst,
                        weight=evtWeight[sel_reg]
                    )
                    output['subjetpt'].fill(
                        jet=msoftdrop,
                        subjetpt=np.minimum(ak.flatten(ak.firsts(fatjet_sortedpnmd.subjets,axis=2).pt[sel_reg]), max(output['subjetpt'].axes[2].edges)-1),
                        systematic=syst,
                        weight=evtWeight[sel_reg]
                    )
                    output['subjet2pt'].fill(
                        jet=msoftdrop,
                        subjet2pt=clip_value(
                            ak.flatten(ak.firsts(fatjet_sortedpnmd.subjets[:,:,1:2],axis=2).pt[sel_reg]),
                            min(output['subjet2pt'].axes[2].edges),
                            max(output['subjet2pt'].axes[2].edges)
                        ),
                        systematic=syst,
                        weight=evtWeight[sel_reg]
                    )
                    output['subjet1_btag'].fill(
                        jet=msoftdrop,
                        subjet1_btag=ak.flatten(ak.firsts(fatjet_sortedpnmd.subjets[:,:,0:1],axis=2).btagDeepB[sel_reg]),
                        systematic=syst,
                        weight=evtWeight[sel_reg]
                    )
                    output['subjet2_btag'].fill(
                        jet=msoftdrop,
                        subjet2_btag=ak.flatten(ak.firsts(fatjet_sortedpnmd.subjets[:,:,1:2],axis=2).btagDeepB[sel_reg]),
                        systematic=syst,
                        weight=evtWeight[sel_reg]
                    )
                    output['subjet1_nBhadrons'].fill(
                        jet=msoftdrop,
                        subjet1_nBhadrons=ak.flatten(ak.firsts(fatjet_sortedpnmd.subjets[:,:,0:1],axis=2).nBHadrons[sel_reg]),
                        systematic=syst,
                        weight=evtWeight[sel_reg]
                    )
                    output['subjet2_nBhadrons'].fill(
                        jet=msoftdrop,
                        subjet2_nBhadrons=ak.flatten(ak.firsts(fatjet_sortedpnmd.subjets[:,:,1:2],axis=2).nBHadrons[sel_reg]),
                        systematic=syst,
                        weight=evtWeight[sel_reg]
                    )
                    output['deltar'].fill(
                        jet=msoftdrop,
                        deltar=ak.flatten(deltar[sel_reg],axis=None),
                        systematic=syst,
                        weight=evtWeight[sel_reg]
                    )
                fatjet_sortedbbvs = FatJet[ak.argsort(FatJet.deepTagMD_bbvsLight, ascending=False)]
                output['hLeadingFatJetDeepTagMD_bbvsLight'].fill(
                    jet=msoftdrop,
                    hLeadingFatJetDeepTagMD_bbvsLight=(ak.flatten(fatjet_sortedbbvs.deepTagMD_bbvsLight[sel_reg][:, 0:1])),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
                )
                output['hLeadingFatJetMass'].fill(
                    jet=msoftdrop,
                    hLeadingFatJetMass=(ak.flatten(fatjet_sortedbbvs.mass[sel_reg][:, 0:1])),
                    systematic=syst,
                    weight=evtWeight[sel_reg]
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
    sInputFiles         = config.pop("inputFiles")
    sOutputFile         = config.pop("outputFile")
    sample_category     = config['sampleCategory']
    process_name        = config.pop('process_name')
    if config['isMC']:
        luminosity          = Luminosities[config['era']][0]
        sample_crossSection = config.pop("crossSection")
        sample_nEvents      = config.pop("nEvents")
        sample_sumEvents    = config.pop("sumEvents") if config["sumEvents"] != -1 else sample_nEvents
        if sample_sumEvents == -1: sample_sumEvents = 1 # Case when sumEvents is not calculated
        lumiScale = calculate_lumiScale(luminosity=luminosity, crossSection=sample_crossSection, sumEvents=sample_sumEvents)
    config['lumiScale'] = lumiScale
    config['isSignal'] = 'SUSY' in process_name
    config['isQCD'] = 'QCD'in process_name
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
            datasetInfo=config
        )
    )
    if 'cutflow' in output.keys():
        print("Cutflow::")
        for msoftdrop in config['msoftdrop']:
            for key in output['cutflow'][msoftdrop].keys():
                print("%10f\t%10d\t%s" % (output['cutflow_weighted'][msoftdrop][key], 
                        output['cutflow'][msoftdrop][key],
                        key)
                  )
    if sOutputFile is not None:
        if not sOutputFile.endswith('.root'): sOutputFile += '.root'
        sDir1 = f'evt/{config["lepton_selection"]}'# % (sample_category)
        with uproot.recreate(sOutputFile) as fOut:
            for key, value in output.items():
                if not isinstance(value, hist.Hist): continue
                for jet in range(0, len(value.axes[0])):
                    jetname = value.axes[0][jet]
                    for syst in  range(0, len(value.axes[1])):
                        systname = value.axes[1][syst]
                        histname = key + '_' + systname if 'central' not in systname else key
                        if 'gen' not in histname:
                            if value.ndim == 3:
                                fOut[f'{sDir1}/{jetname}/{sample_category}/{histname}'] = value[jet,syst,:]
                            else:
                                fOut[f'{sDir1}/{jetname}/{sample_category}/{histname}'] = value[jet,syst,:,:]
                        else:
                            fOut[f'{sDir1}/{jetname}/{sample_category}/gen/{histname}'] = value[jet,syst,:]
    current_memory, peak_memory = tracemalloc.get_traced_memory() # https://medium.com/survata-engineering-blog/monitoring-memory-usage-of-a-running-python-program-49f027e3d1ba
    print(f"\n\nMemory usage:: current {current_memory / 10**6}MB;  peak {peak_memory / 10**6}MB")

    endTime = time.time()
    totalTime = endTime - startTime
    totalTime_hr  = int(totalTime/60/60)
    totalTime_min = totalTime - float(totalTime_hr * 60)
    totalTime_min = int(totalTime_min/60)
    totalTime_sec = totalTime - float(totalTime_hr * 60*60) - float(totalTime_min * 60)
    #for f in sInputFiles:
    #os.system(f'rm {f}')
    print(f"Total run time: {totalTime_hr}h {totalTime_min}m {totalTime_sec}s = {totalTime}sec ")
