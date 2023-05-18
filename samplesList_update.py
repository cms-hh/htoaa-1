import uproot3 as uproot
import numpy as np
import awkward as ak
import sys
import json

encoding = 'utf-8'
old_file  = 'Samples_2018UL_v0.json'
with open(old_file) as fSamplesInfo:
        samplesInfo = json.load(fSamplesInfo) # Samples_Era.json

def get_histvalue(proc):
    with uproot.open('hadd_countgenWeight.root') as f:
        for key in f['evt'].keys():
            proc_name = str(f['evt'][key].name, encoding)
            if proc_name != proc: continue
            hist = f[f'evt/{proc_name}/hCutFlow_central']
            return hist.values[1]

def check_sample(proc, samples_sum):
        for sample_sum in samples_sum:
                if proc in sample_sum[0]: return True
        return False

samples_sum = []
for idx1, sample_info in enumerate(samplesInfo):
        samplesInfo[sample_info]['sumEvents'] = get_histvalue(samplesInfo[sample_info]['process_name'])

with open('Samples_2018UL.json', "w") as fSampleInfo:
        json.dump(samplesInfo, fSampleInfo, indent=4)
