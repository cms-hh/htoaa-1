import uproot as uproot
import numpy as np
import awkward as ak
import hist
import argparse
import os
import sys
import json
import time
import tracemalloc
from htoaa_CommonTools import GetDictFromJsonFile, get_lf
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema, BaseSchema
from coffea import processor
from coffea.nanoevents import schemas


nEventToReadInBatch = 0.5*10**6 # 2500000 #  1000 # 2500000
nEventsToAnalyze =  -1 #-1 # 1000 # 100000 # -1

class PileupProcessor(processor.ProcessorABC):

    def __init__(self, datasetInfo={}):

        self.pu_data=uproot.open('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/UltraLegacy/PileupHistogram-goldenJSON-13tev-2018-69200ub-99bins.root')['pileup'].to_hist()

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = {}
        output['pu_mc'] = (hist.Hist.new
                           .Reg(len(self.pu_data.values()),
                                self.pu_data.axes.edges[0][0],
                                self.pu_data.axes.edges[0][len(self.pu_data.axes.edges[0])-1],
                                name='pu_mc', label='pu_mc')
                           .Weight()
                       )
        output['pu_mc'].fill(events.Pileup.nTrueInt)
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

    sOutputFile         = config["outputFile"]
    sample_category     = config['sampleCategory']
    process_name        = config['process_name']
    assert(config["isMC"])
    era                 = config['era']
    sInputFiles = get_lf(config["inputFiles"], process_name)
    startTime = time.time()
    tracemalloc.start()
    chunksize = nEventToReadInBatch
    maxchunks = None if nEventsToAnalyze == -1 else int(nEventsToAnalyze/nEventToReadInBatch)
    nWorkers  = 4 if nEventsToAnalyze == -1 else 1
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
        processor_instance=PileupProcessor(
            datasetInfo=config
            )
    )
    if sOutputFile is not None:
        if not sOutputFile.endswith('.root'): sOutputFile += '.root'
        sDir1 = 'evt/%s' % (process_name) #
        
        with uproot.recreate(sOutputFile) as fOut:
            for key, value in output.items():
                if not isinstance(value, hist.Hist): continue
                #print(value.shape)
                fOut[f'{sDir1}/{process_name}'] = value
                
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
'''
pu_data=uproot.open('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/UltraLegacy/PileupHistogram-goldenJSON-13tev-2018-69200ub-99bins.root')['pileup'].to_hist()

events = NanoEventsFactory.from_root(
    "/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/PNet_v1_2023_10_06/SUSY\
_WH_WToAll_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8/skims/Hto4b_0p8/PNet_v1_skim_Hto4b_0p8_WH_WToAll_M-12.root",
    schemaclass=NanoAODSchema.v6
).events()
pu_mc = (hist.Hist.new
         .Variable(pu_data.axes.edges[0], name='pu_mc', label='pu_mc')
         .Weight()
     ).fill(events.Pileup.nTrueInt)
pu_w=np.divide(pu_data.values()/pu_data.values().sum(), pu_mc.values()/pu_mc.values().sum(), where=pu_mc.values()>0, out=np.ones_like(pu_data.values()))
while True:
    x= np.minimum(x,3)
    print(x)
    print('p: ', (np.sum(x)-np.sum(y))/np.sum(y))
    sf = 1/(1+(np.sum(x)-np.sum(y))/np.sum(y))
    x *= sf
    print('sf: ', sf)
    if sf <= (1+0.0025):
        break
'''
