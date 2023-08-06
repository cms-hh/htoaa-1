#htoaa analysis main code

import os
import sys
import subprocess
import json
import glob
from collections import OrderedDict as OD
import time
import tracemalloc
import math
import numpy as np
from copy import copy, deepcopy
#import uproot
import uproot3 as uproot

'''
Count sum of positive and negative gen-weight events in sample

References:
  * Coffea framework used for TTGamma analysis: https://github.com/nsmith-/TTGamma_LongExercise/blob/FullAnalysis/ttgamma/processor.py
  * Coffea installation: /home/siddhesh/anaconda3/envs/ana_htoaa/lib/python3.10/site-packages/coffea
'''

#import coffea.processor as processor
from coffea import processor, util
from coffea.nanoevents import schemas
from coffea.nanoevents.methods import nanoaod, vector
from coffea.analysis_tools import PackedSelection, Weights
#import hist
from coffea import hist
import awkward as ak
import uproot
#from dask.distributed import Client

from htoaa_Settings import *
from htoaa_CommonTools import (
    GetDictFromJsonFile, setXRootDRedirector,
    xrdcpFile
)

printLevel = 0
nEventToReadInBatch =  0.5*10**6 # 2500000 #  1000 # 2500000
nEventsToAnalyze = -1 # 1000 # 100000 # -1


sWeighted = "Wtd: "

class HToAATo4bProcessor(processor.ProcessorABC):
    def __init__(self, datasetInfo={}):

        ak.behavior.update(nanoaod.behavior)
        self.datasetInfo = datasetInfo
        dataset_axis    = hist.Cat("dataset", "Dataset")
        systematic_axis = hist.Cat("systematic", "Systematic Uncertatinty")
        cutFlow_axis  = hist.Bin("CutFlow",   r"Cuts", 21, -0.5, 20.5)
        sXaxis      = 'xAxis'
        sXaxisLabel = 'xAxisLabel'
        sYaxis      = 'yAxis'
        sYaxisLabel = 'yAxisLabel'
        histos = OD([
            ('hCutFlow', {sXaxis: cutFlow_axis,    sXaxisLabel: 'Cuts'}),
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
                        "Counts",
                        dataset_axis,
                        hXaxis, #nObject_axis,
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
                        hXaxis, #nObject_axis,
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

        if nEventsToAnalyze != -1:
            print(f"\n (run:ls:event): {ak.zip([events.run, events.luminosityBlock, events.event])}")
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

        ones_list  = np.ones(len(events))
        selection = PackedSelection()

        sel_posGenWgt = None
        sel_negGenWgt = None        
        
        if self.datasetInfo[dataset]["isMC"]:
            selection.add("posGenWgt", (events.genWeight > 0))
            selection.add("negGenWgt", (events.genWeight < 0))

            sel_posGenWgt           = selection.all("posGenWgt")
            sel_negGenWgt           = selection.all("negGenWgt")
        
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
        if self.datasetInfo[dataset]["isMC"]:
            weights_gen.add(
                "genWeight",
                weight=events.genWeight#np.copysign(np.ones(len(events)), events.genWeight)
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
            else:
                systList = [shift_syst]
        else:
            systList = ["noweight"]

        output['cutflow']['all events'] += len(events)
        output['cutflow'][sWeighted+'all events'] += weights.weight().sum() 
        for n in selection.names:
            sel_i = selection.all(n)
            output['cutflow'][n] += sel_i.sum()
            output['cutflow'][sWeighted+n] += weights.weight()[sel_i].sum()

        for syst in systList:

            # find the event weight to be used when filling the histograms
            weightSyst = syst
            
            # in the case of 'central', or the jet energy systematics, no weight systematic variation is used (weightSyst=None)
            if syst in ["central", "JERUp", "JERDown", "JESUp", "JESDown"]:
                weightSyst = None
            if syst == "noweight":
                evtWeight = np.ones(len(events))
            else:
                evtWeight     = weights.weight(weightSyst)
                evtWeight_gen = weights_gen.weight(weightSyst)
            # all events
            iBin = 0
            output['hCutFlow'].fill(
                dataset=dataset,
                CutFlow=(ones_list * iBin),
                systematic=syst
            )
            # MC ----------------------------------------------
            iBin = 1
            if self.datasetInfo[dataset]["isMC"]:
                output['hCutFlow'].fill(
                    dataset=dataset,
                    CutFlow=(ones_list * iBin),
                    systematic=syst,
                    weight=evtWeight_gen
                )
            else:
                output['hCutFlow'].fill(
                    dataset=dataset,
                    CutFlow=(ones_list * iBin),
                    systematic=syst,
                    weight=weights
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
    #branchesToRead = htoaa_nanoAODBranchesToRead
    #print("branchesToRead: {}".format(branchesToRead))

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
    startTime = time.time()
    tracemalloc.start()
    chunksize = nEventToReadInBatch
    maxchunks = None if nEventsToAnalyze == -1 else int(nEventsToAnalyze/nEventToReadInBatch)
    nWorkers  = 4 if nEventsToAnalyze == -1 else 1
    print(f"nEventsToAnalyze: {nEventsToAnalyze},  nEventToReadInBatch: {nEventToReadInBatch}, chunksize: {chunksize},  maxchunks: {maxchunks},  nWorkers: {nWorkers}")
    run = processor.Runner(
        #executor=executor,
        executor=processor.FuturesExecutor(workers=nWorkers),
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
                sample_category: {"isMC": isMC}
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
        sDir1 = 'evt/%s' % (process_name) #
        
        with uproot.recreate(sOutputFile) as fOut:
            for key, value in output.items():
                if not isinstance(value, hist.Hist): continue
                for _dataset in value.axis('dataset').identifiers():
                    for _syst in value.axis('systematic').identifiers():
                        h1 = value.integrate('dataset',_dataset).integrate('systematic',_syst).to_hist()
                        fOut['%s/%s_%s' % (sDir1, key, _syst)] = h1
                
        print("Wrote to sOutputFile {}".format(sOutputFile))

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
