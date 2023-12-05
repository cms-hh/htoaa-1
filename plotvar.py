import matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
matplotlib.use('agg')
import uproot3 as uproot
import numpy as np
from math import sqrt
import awkward as ak
import hist
import argparse
import os

parser = argparse.ArgumentParser(description='htoaa analysis wrapper')
parser.add_argument('-p', dest='input_path', type=str, required=True)
parser.add_argument('-o', dest='output_file', type=str, required=True)
parser.add_argument('-v', dest='variables', nargs='+', default=[])
parser.add_argument('-b', dest='backgrounds', nargs='+', default=[])
parser.add_argument('-s', dest='signals', nargs='+', default=[])
args=parser.parse_args()

p = args.input_path
output_path = args.output_file
variables = args.variables
backgrounds = args.backgrounds
signals = args.signals

f=uproot.open(p)
if not os.path.exists(output_path): os.mkdir(output_path)

for var in variables:

    hists = {}
    for bkg in backgrounds:
        hists[bkg] = f[f'evt/tight_lep/msoftdrop_GE_90/{bkg}'][var]

    hists = sorted(hists.items(), key=lambda h: np.sum(h[1].values), reverse=False)
    sorted_hists = {}
    for h in hists:
        sorted_hists[h[0]] = h[1]

    sig_hists = {}
    for sig in signals:
        sig_hists[sig] = f[f'evt/tight_lep/msoftdrop_GE_90/{sig}'][var]
    
    #fig, axis = plt.subplots(ncols=1, nrows=2, figsize=(8,10), sharex='col', gridspec_kw={'height_ratios': [3, 1], 'hspace': 0})
    fig, axis = plt.subplots(figsize=(10,10), sharex='col')
    max_bkg = np.max([np.max(sorted_hists[k].values) for k in sorted_hists.keys()])
    max_sig = np.max([np.max(sig_hists[k].values) for k in sig_hists.keys()])
    max_v = np.max([max_bkg, max_sig])
    #axis0=axis[0]
    axis.semilogy()
    hep.histplot([sorted_hists[v].values for v in sorted_hists.keys()],
                 bins=sorted_hists['TT'].edges,
                 ax=axis, label=[v for v in sorted_hists.keys()], histtype='step')
    hep.histplot([sig_hists[v].values for v in sig_hists],
                 bins=sorted_hists['TT'].edges,
                 ax=axis, label=[v for v in sig_hists.keys()])
    axis.legend(fontsize=12, loc='best', ncol=2, bbox_to_anchor=(-0.1, 0.65, 1.1, 0.36))
    print(max_v)
    axis.set_ylim(0.001, 10*max_v)
    #plt.tight_layout()
    plt.savefig(os.path.join(output_path, f'{var}.png'))
    print(output_path)
    #plt.savefig(output_path)
