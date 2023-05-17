import uproot3 as uproot
import numpy as np
import awkward as ak
import sys
f=uproot.open('hadd.root')
bkgs = ['WJets', 'TT']
sigs = ['SUSY_WH_WToAll_HToAATo4B_Pt150_M-12', 'SUSY_WH_WToAll_HToAATo4B_Pt150_M-15', 'SUSY_WH_WToAll_HToAATo4B_Pt150_M-20', 'SUSY_WH_WToAll_HToAATo4B_Pt150_M-25', 'SUSY_WH_WToAll_HToAATo4B_Pt150_M-30', 'SUSY_WH_WToAll_HToAATo4B_Pt150_M-35', 'SUSY_WH_WToAll_HToAATo4B_Pt150_M-40', 'SUSY_WH_WToAll_HToAATo4B_Pt150_M-45', 'SUSY_WH_WToAll_HToAATo4B_Pt150_M-50', 'SUSY_WH_WToAll_HToAATo4B_Pt150_M-55', 'SUSY_WH_WToAll_HToAATo4B_Pt150_M-60']
vars = ['hLeadingFatJetParticleNetMD_Xbb_central', 'hLeadingFatJetDeepTagMD_bbvsLight_central']
sighist = f['evt']['SUSY_WH_WToAll_HToAATo4B_Pt150_M-12'][vars[0]]
bincenters = (sighist.edges[:-1] + sighist.edges[1:]) / 2
x = np.linspace(0,1,100)
for var in vars:
    for sig in sigs:
        sOb = []
        for ival, val in enumerate(bincenters):
            sighist = f['evt'][sig][var]
            sigint = ak.sum(np.where(bincenters>=val, sighist.values, ak.zeros_like(sighist.values)))
            for idx, bkg in enumerate(bkgs):
                bkghist = f['evt'][bkg][var]
                if idx == 0:
                    bkgint = ak.sum(np.where(bincenters>=val, bkghist.values, ak.zeros_like(bkghist.values)))
                else:
                    bkgint += ak.sum(np.where(bincenters>=val, bkghist.values, ak.zeros_like(bkghist.values)))
            bkgint = np.sqrt(bkgint)
            sOb.append(sigint/np.sqrt(bkgint))
        sOb = np.array(sOb)
        sOb = np.nan_to_num(sOb, nan=0).tolist()
        print(var)
        print(sig)
        print(bincenters[sOb.index(max(sOb))])
        if sOb.index(max(sOb)) != len(bincenters)-1:
            print(f"{bincenters[sOb.index(max(sOb))]:0,.3f}, \t {sOb[sOb.index(max(sOb))]:0,.2f}")
        else:
            print(f"{bincenters[sOb.index(max(sOb))]:0,.3f}, \t {sOb[sOb.index(max(sOb))]:0,.2f}")
