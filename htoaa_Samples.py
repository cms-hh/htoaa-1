#import OrderedDict as OD
from collections import OrderedDict as OD

kData = "data_obs" # dict key for Datasets

Samples2018 = OD([

#    (kData, [
#        "/JetHT/Run2018A-UL2018_MiniAODv2_NanoAODv9-v2/NANOAOD",
#        "/JetHT/Run2018B-UL2018_MiniAODv2_NanoAODv9-v1/NANOAOD",
#        "/JetHT/Run2018C-UL2018_MiniAODv2_NanoAODv9-v1/NANOAOD",
#        "/JetHT/Run2018D-UL2018_MiniAODv2_NanoAODv9-v2/NANOAOD"
#    ]),    
    (kData, [
        "JetHT_Run2018A",
        "JetHT_Run2018B",
        "JetHT_Run2018C",
        "JetHT_Run2018D"
    ]),
    
    ("QCD_bEnrich",[
        "QCD_bEnriched_HT100to200_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_bEnriched_HT200to300_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_bEnriched_HT300to500_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_bEnriched_HT500to700_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_bEnriched_HT700to1000_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_bEnriched_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_bEnriched_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_bEnriched_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8"
    ]),

    ("QCD_bGen",[
        "QCD_HT100to200_BGenFilter_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_HT200to300_BGenFilter_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_HT300to500_BGenFilter_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_HT500to700_BGenFilter_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_HT700to1000_BGenFilter_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_HT1000to1500_BGenFilter_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_HT1500to2000_BGenFilter_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_HT2000toInf_BGenFilter_TuneCP5_13TeV-madgraph-pythia8"
    ]),

    ("QCDIncl",[
        "QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8",
        "QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8",
        "QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8",
        "QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8",
        "QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8",
        "QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8",
        "QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8",
        "QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8",
        "QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8"
    ]),

#    ("QCDIncl_PSWeight", [
#        "QCD_HT50to100_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
#        "QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
#        "QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
#        "QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
#        "QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
#        "QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
#        "QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
#        "QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
#        "QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
#     ]),

    ("ZJets", [
        "ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8",
        "ZJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8",
        "ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8",
        "ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8"
    ]),
    
    ("WJets", [
        "WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8"
    ]),
    

    ("SUSY_GluGluH_01J_HToAATo4B", [
        #"SUSY_GluGluH_01J_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8",
        "SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8"
    ]),
    
    ("SUSY_WH_WToAll_HToAATo4B", [
        "SUSY_WH_WToAll_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8",
        "SUSY_WH_WToAll_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8",
        "SUSY_WH_WToAll_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8",
        "SUSY_WH_WToAll_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8",
        "SUSY_WH_WToAll_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8",
        "SUSY_WH_WToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8",
        "SUSY_WH_WToAll_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8",
        "SUSY_WH_WToAll_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8",
        "SUSY_WH_WToAll_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8",
        "SUSY_WH_WToAll_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8",
        "SUSY_WH_WToAll_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8"
    ]),

    ('WJetsToLNu_HT', [
        "WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8"
    ]),

    ("TT", [
        #"TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8",
        "TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8",
        "TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8"
    ])
])
