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


printLevel = 0

sXS= "xs"
sNameSp = "nameSp"
list_datasetAndXs_2018 = OD([
    
    ## QCD_bEnriched_HT*
    ("/QCD_bEnriched_HT100to200_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:1117000.0}),
    ("/QCD_bEnriched_HT200to300_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:80430.0}),
    ("/QCD_bEnriched_HT300to500_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:16620.0}),
    ("/QCD_bEnriched_HT500to700_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:1487.0}),
    ("/QCD_bEnriched_HT700to1000_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:296.5}),
    ("/QCD_bEnriched_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:46.61}),
    ("/QCD_bEnriched_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:3.72}),
    ("/QCD_bEnriched_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:0.6462}),
    #dasgoclient --query="dataset=/QCD_bEnriched_HT*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/QCD_HT100to200_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:1275000.0}),
    ("/QCD_HT200to300_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:111700.0}),
    ("/QCD_HT300to500_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:27960}),
    ("/QCD_HT500to700_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:3078}),
    ("/QCD_HT700to1000_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:721.8}),
    ("/QCD_HT1000to1500_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:138.2}),
    ("/QCD_HT1500to2000_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:13.61}),
    ("/QCD_HT2000toInf_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:2.92}),
    ## QCD_HT*_BGenFilter
    # dasgoclient --query="dataset=/QCD_HT*_BGenFilter*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ## QCD_HT* TuneCP5 madgraphMLM
    # dasgoclient --query="dataset=/QCD_HT*TuneCP5_13TeV-madgraphMLM*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:185300000}),
    ("/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:23590000}),
    ("/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:1551000}),
    ("/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:323400}),
    ("/QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:30140}),
    ("/QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:6344}),
    ("/QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:1092}),
    ("/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:99.76}),
    ("/QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:20.35}),
    ## QCD_HT* PSWeights
    # dasgoclient --query="dataset=/QCD_HT*PSWeight*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/QCD_HT50to100_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:187300000}),
    ("/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:23590000}),
    ("/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:1555000}),
    ("/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:324500}),
    ("/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:30310}),
    ("/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:6444}),
    ("/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:1127}),
    ("/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:109.8}),
    ("/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:21.98}),
    ## TTbar Jets
    # dasgoclient --query="dataset=/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    #("/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 831.76}),
    ("/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:109.6}),
    ("/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:108.7}),
    ("/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:54.23}),
    


    ## ZJets
    # dasgoclient --query="dataset=/ZJetsToQQ*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:1012}),
    ("/ZJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:114.2}),
    ("/ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:25.34}),
    ("/ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:190.2 }),
    ## WJets
    # dasgoclient --query="dataset=/WJetsToQQ*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:276.5 }),
    ("/WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:59.25 }),
    ("/WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:28.75 }),
    ("/W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:10167.853 }),
    ("/W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:3199.4510 }),
    ("/W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:941.15570 }),
    ("/W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:439.084 }),
    ("/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:1440}),
    ("/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 382.1}), #https://cms-gen-dev.cern.ch/xsdb/?columns=67108863&currentPage=0&pageSize=10&searchQuery=DAS%3DTTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8 3 xs
    ("/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:61526.7}),
    ("/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 51.54}), #https://cms-gen-dev.cern.ch/xsdb/?columns=67108863&currentPage=0&pageSize=10&searchQuery=DAS%3DWJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8
    ("/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM", {sXS: 51.54}),
    ("/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 12.49}),
    ("/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM", {sXS:12.49}),
    ("/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 5.619}),
    ("/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM", {sXS:5.619}),
    ("/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1.321}),
    ("/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM", {sXS:1.321}),
    ("/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS: 0.02992}),
    ("/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1431}),
    
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

        datasetName_parts            = datasetName.split('/')[1:]
        sampleName                   = datasetName_parts[0]
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
