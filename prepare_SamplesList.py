import subprocess
import json
from collections import OrderedDict as OD
from copy import deepcopy
import argparse
import os
import glob
import uproot
import sys
from htoaa_Settings import *


printLevel = 0

sXS= "xs"
sNameSp = "nameSp"
list_datasetAndXs_2018 = OD([
    
    ## QCD_bEnriched_HT*
    # dasgoclient --query="dataset=/QCD_bEnriched_HT*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"    
    ## QCD_HT*_BGenFilter
    # dasgoclient --query="dataset=/QCD_HT*_BGenFilter*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ## QCD_HT* TuneCP5 madgraphMLM
    # dasgoclient --query="dataset=/QCD_HT*TuneCP5_13TeV-madgraphMLM*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ## QCD_HT* PSWeights
    # dasgoclient --query="dataset=/QCD_HT*PSWeight*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ## TTbar Jets
    # dasgoclient --query="dataset=/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 831.76}),
    ("/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9/NANOAODSIM", {sXS:108.7}),
    ("/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9/NANOAODSIM", {sXS:109.6}),
    


    ## ZJets
    # dasgoclient --query="dataset=/ZJetsToQQ*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"    
    ## WJets
    # dasgoclient --query="dataset=/WJetsToQQ*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9/NANOAODSIM", {sXS: 407.9}), #https://cms-gen-dev.cern.ch/xsdb/?columns=67108863&currentPage=0&pageSize=10&searchQuery=DAS%3DTTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8 3 xs
    ("/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9/NANOAODSIM", {sXS: 57.48}), #https://cms-gen-dev.cern.ch/xsdb/?columns=67108863&currentPage=0&pageSize=10&searchQuery=DAS%3DWJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8
    ("/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9/NANOAODSIM", {sXS: 12.87}),
    ("/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9/NANOAODSIM", {sXS: 5.366}),
    ("/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9/NANOAODSIM", {sXS: 1.074}),
    ("/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9/NANOAODSIM", {sXS: 0.008001}),
    ("/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9/NANOAODSIM", {sXS: 1395}),
    
    ## SUSY_GluGluH_01J_HToAATo4B_M-*   and   SUSY_GluGluH_01J_HToAATo4B_Pt150_M-*
    # dasgoclient --query="dataset=/SUSY*GluGluH*HToAATo4B*M*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ## JetHT
    # dasgoclient --query="dataset=/JetHT/*2018*UL*MiniAODv2_NanoAODv9-*/NANOAOD"
    # XS (cross-section) does not matter for data sample
    #WH
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1.}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1}),
])


sDataset = "dataset"
sNanoAOD_nFiles = "nanoAOD_nFiles"
sNanoAOD = "nanoAOD"
sCross_section = "cross_section"
sNEvents = "nEvents"
sSumEvents = "sumEvents"
sc = "sample_category"
pn = "process_name"
sampleDetail_dict_template = OD([
    (sDataset, ""),
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
    for idx, file in enumerate(files):
        with uproot.open(file) as f:
            nEventsTotal += len(f['Events']['event'].array())
    return nEventsTotal, nFiles, files

def getprocess_name(name):
    ret = name.split('_Tune')[0]
    return ret
    if name.startswith('TT'):
        return 'TT'
    elif name.startswith('WJets'):
        return 'WJets'
    elif name.startswith('SUSY'):
        return name

def getsample_category(name):
    ret = name.split('_Tune')[0].split('_exit')[0]
    if 'PSweights' in name:
        ret += '_PSweights'
    #if 

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
    args=parser.parse_args()

    era           = args.era
    print(f"era: {era}")


    list_datasetAndXs = None
    sFileSamplesInfo_toUse = None
    if era == Era_2018:
        list_datasetAndXs = list_datasetAndXs_2018
        sFileSamplesInfo_toUse = sFileSamplesInfo[era]

    sFileSamplesInfo_toUse = sFileSamplesInfo_toUse.replace('.json', '_v0.json')

        
    samples_details = OD()
    for datasetName, datasetDetails in list_datasetAndXs.items():
        sampleDetails_dict = deepcopy(sampleDetail_dict_template)
        #print(f"sampleDetails_dict_0: {sampleDetails_dict}")

        sampleDetails_dict[sDataset] = datasetName
        datasetName_parts            = datasetName.split('/')[1:]
        print('datasetName_parts: ', datasetName_parts)
        sampleName                   = datasetName_parts[0]
        
        if datasetName_parts[-1] == 'NANOAODSIM':
            # for MC sample
            sampleDetails_dict[sCross_section] = datasetDetails[sXS]
        else:
            # for data sample
            sampleName_part2 = (datasetName_parts[1]).split('-')[0] # 'Run2018A-UL2018_MiniAODv2_NanoAODv9-v2'
            sampleName = '%s_%s' % (sampleName, sampleName_part2)  # JetHT_Run2018A
            del sampleDetails_dict[sCross_section]
            del sampleDetails_dict[sSumEvents]

        
        
        nEventsTotal, nFiles, files = getDatasetFilesFromeos(datasetName_parts)#getDatasetFiles(datasetName)
        sampleDetails_dict[sNEvents] = nEventsTotal
        sampleDetails_dict[sNanoAOD_nFiles] = nFiles
        sampleDetails_dict[pn] = getprocess_name(datasetName_parts[0])
        sampleDetails_dict[sc] = getsample_category(datasetName_parts[0])
        sampleDetails_dict[sNanoAOD] = files
        samples_details[sampleName] = sampleDetails_dict

    if printLevel >= 0:
        print("\n\nsamples:: \n",json.dumps(samples_details, indent=4))

    with open(sFileSamplesInfo_toUse, "w") as fSampleInfo:
        json.dump(samples_details, fSampleInfo, indent=4)

        print(f"\n\nNew sample list wrote to {sFileSamplesInfo_toUse}")
