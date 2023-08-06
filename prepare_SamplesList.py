import subprocess
import json
from collections import OrderedDict as OD
from copy import deepcopy
import argparse
import os
import glob
import uproot
import numpy as np
import sys
import re
from htoaa_Settings import *
import importlib

printLevel = 0

sXS= "xs"
sNameSp = "nameSp"
sNanoAOD_nFiles = "nanoAOD_nFiles"
sNanoAOD = "nanoAOD"
sCross_section = "cross_section"
sNEvents = "nEvents"
sSumEvents = "sumEvents"
sc = "sample_category"
pn = "process_name"
sampleDetail_dict_template = OD([
    (sCross_section, -1.),
    (sNEvents, -1),
    (sSumEvents, -1),
    (sNanoAOD, []),
    (sc, ''),
    (sNanoAOD_nFiles, -1)
])
#/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9/112F42A6-4315-D547-9161-DB39ECC06CB2.root
def getDatasetFilesFromeos(datasetName_parts):
    prefix = '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC'
    campaign = datasetName_parts[1].split('-')[0]
    path = os.path.join(prefix, datasetName_parts[0], campaign)
    files = glob.glob(f'{path}/*root')
    nFiles = len(files)
    nEventsTotal = 0
    sumEvents = 0
    for idx, file in enumerate(files):
        with uproot.open(file) as f:
            nevents = np.array(f['Events']['event'].array())
            nEventsTotal += len(nevents)
            nevents.fill(1)
            genevents = np.array(f['Events']['genWeight'].array())
            nevents = np.copysign(nevents, genevents)
            sumEvents += np.sum(nevents)
    return sumEvents, nEventsTotal, nFiles, files

def getprocess_name(datasetName_parts):
    campaign = datasetName_parts[1]#.split('-')[0]
    name = datasetName_parts[0]
    ret = name.split('_Tune')[0]
    if 'PSWeights' in name:
        ret += '_PSWeights'
    if 'ext' in campaign:
        ext, v = campaign.split('-')[-2], campaign.split('-')[-1]
        ext = ext.split('_')[-1]
        ret += '_'+ext
        ret += '-'+v
    return ret

def getsample_category(name):
    name = name.split('_Tune')[0]
    if name.startswith('TT'):
        return 'TT'
    elif name.startswith('WJet') or re.findall('^W[0-4]Jets', name):
        return 'WJets'
    elif name.startswith('SUSY'):
        return name
    elif name.startswith('QCD'):
        return 'QCD'
    elif name.startswith('ZZ') or name.startswith('ZJets'):
        return 'Z'
    else:
        print(name)
        assert(0)

def getDatasetFiles(dataset):
    print('dataset: ', dataset)
    cmd1 = ['bash','-c', 'dasgoclient --query="file dataset=%s" --format=json'%(dataset)]
    if printLevel >= 10:
        print(f"cmd1: {cmd1}")
    output = subprocess.check_output(cmd1) # output in bytes
    output = output.decode("utf-8") # convert output in bytes to string
    output = json.loads(output) # convert output in 'string' to 'dict'
    nFiles = output['nresults']
    files  = []
    nEventsTotal = 0

    if nFiles != len(output['data']):
        print(f"nFiles != len(output['data']... something is wrong.. \t\t **** ERROR ****")
        exit(0)
    for iFile in range(nFiles):
        if len(output['data'][iFile]['file']) != 1:
            print(f"len(output['data'][iFile]['file']) != 1: Assumption of a single entry list 'output['data'][iFile]['file']' seems wrong... need to follow up.. \t\t **** ERROR **** ")
            exit(0)
        file_LFN = output['data'][iFile]['file'][0]['name']
        nEvents  = output['data'][iFile]['file'][0]['nevents']
        if printLevel >= 5:
            print(f"file_LFN: {file_LFN}, nEvents ({type(nEvents)}): {nEvents}, nEventsTotal: {nEventsTotal}  {output['data'][iFile]['file'][0]}")
            
        files.append( file_LFN )
        nEventsTotal += nEvents

    if printLevel >= 3:
        print(f"\ndataset: {dataset}, nEventsTotal: {nEventsTotal}, nFiles: {nFiles}, files: {files}")
    return nEventsTotal, nFiles, files
    

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='htoaa analysis wrapper')
    parser.add_argument('-era', dest='era', type=str, default=Era_2018, choices=[Era_2016, Era_2017, Era_2018], required=False)

    args = parser.parse_args()
    era = args.era

    list_datasetAndXs = None
    sFileSamplesInfo_toUse = None
    if era == Era_2018:
        mod = importlib.import_module('samplesList_2018')
        list_datasetAndXs = getattr(mod,'list_datasetAndXs')
        #list_datasetAndXs = list_datasetAndXs_2018
        sFileSamplesInfo_toUse = sFileSamplesInfo[era]

    sFileSamplesInfo_toUse = sFileSamplesInfo_toUse.replace('.json', '_v0.json')

    samples_details = OD()
    for datasetName, datasetDetails in list_datasetAndXs.items():
        sampleDetails_dict = deepcopy(sampleDetail_dict_template)
        datasetName_parts = datasetName.split('/')[1:]
        sampleName = datasetName_parts[0]
        if sampleName in samples_details.keys():
            print(f'same sample {sampleName} is used')
            assert(0)
        if datasetName_parts[-1] == 'NANOAODSIM':
            # for MC sample
            sampleDetails_dict[sCross_section] = datasetDetails[sXS]
        else:
            # for data sample
            sampleName_part2 = (datasetName_parts[1]).split('-')[0] # 'Run2018A-UL2018_MiniAODv2_NanoAODv9-v2'
            sampleName = '%s_%s' % (sampleName, sampleName_part2)  # JetHT_Run2018A
            del sampleDetails_dict[sCross_section]
            del sampleDetails_dict[sSumEvents]

        sumEvents = -1
        if 0:
            sumEvents, nEventsTotal, nFiles, files = getDatasetFilesFromeos(datasetName_parts)
        else:
            nEventsTotal, nFiles, files = getDatasetFiles(datasetName)
        sampleDetails_dict[sNEvents] = nEventsTotal
        sampleDetails_dict[sSumEvents] = sumEvents
        sampleDetails_dict[sNanoAOD_nFiles] = nFiles
        sampleDetails_dict[pn] = getprocess_name(datasetName_parts)
        sampleDetails_dict[sc] = getsample_category(datasetName_parts[0])
        sampleDetails_dict[sNanoAOD] = files
        samples_details[datasetName] = sampleDetails_dict

    if printLevel >= 0:
        print("\n\nsamples:: \n",json.dumps(samples_details, indent=4))

    with open(sFileSamplesInfo_toUse, "w") as fSampleInfo:
        json.dump(samples_details, fSampleInfo, indent=4)
        print(f"\n\nNew sample list wrote to {sFileSamplesInfo_toUse}")
