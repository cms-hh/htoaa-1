import uproot3 as uproot
import sys
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
import mplhep as hep
import numpy as np

file_='hadd.root'
encoding = 'utf-8'
for var in ['hLeadingFatJetParticleNetMD_Xbb_central', 'hLeadingFatJetDeepTagMD_bbvsLight_central']:
    d = uproot.open('hadd.root')['evt']
    foundbkg = False
    for key in d.keys():
        name = str(d[key].name, encoding)
        if name == 'TTJets': continue
        if 'SUSY' not in name:
            if not foundbkg:
                bkg = d[f'{name}/{var}'].values
                err2 = d[f'{name}/{var}'].variances
                foundbkg = True
            else:
                bkg = np.sum((bkg, d[f'{name}/{var}'].values), axis=0)
                err2 = np.sum((err2, d[f'{name}/{var}'].variances), axis=0)
    for m in [12,15,20,25,30,35,40,45,50,55,60]:
        for key in d.keys():
            name = str(d[key].name, encoding)
            if 'SUSY' not in name: continue
            if str(m) not in name:continue
            hist = d[f'{name}/{var}']
            plt.style.use([hep.style.CMS])
            fig, axis = plt.subplots()
            hep.histplot(hist.values, bins=hist.edges, ax=axis, yerr=np.sqrt(hist.variances), density=True, label=f'mass_{str(m)}')
            hep.histplot(bkg, bins=hist.edges, ax=axis, yerr=np.sqrt(err2), density=True, label=f'bkg')
            axis.set_yscale('log', base=10)
            plt.legend()
            plt.tight_layout()
            plt.savefig(f"{var}_mass_{str(m)}.pdf")
            plt.close('all')