from ROOT import TH1F, TH1D, TH1, TFile
import os
import re
import argparse
import glob
import sys
parser = argparse.ArgumentParser(description='addBackgrounds')
parser.add_argument('-p', dest='input_path', type=str, required=True)
parser.add_argument('-o', dest='output_file')#, type=str, required=True)
parser.add_argument('-s', dest='samples_list', type=str, nargs='+', required=True)
parser.add_argument('-syst', dest='systematic_list', type=str, nargs='+', default=[])
args=parser.parse_args()
p = args.input_path
output_file = args.output_file
fo = TFile(output_file, 'recreate')
sl = args.samples_list
syl = args.systematic_list
fs = TFile(p)

for msoftdrop in ['msoftdrop_L_90', 'msoftdrop_GE_90']:
  fo.mkdir(f'evt/tight_lep/{msoftdrop}/nonQCD')
  fo.mkdir(f'evt/tight_lep/{msoftdrop}/data_QCD')
  d = fs.Get('evt')
  sd = d.Get(f'tight_lep/{msoftdrop}/{sl[0]}')
  print(sd, d)
  data_obs = fs.Get(f'evt/tight_lep/{msoftdrop}/data_obs')
  histnames_list = [h.GetName() for h in sd.GetListOfKeys()]
  histnames_list_wosyst = []
  for histname in histnames_list:
    contain_syst = False
    for syst in syl:
      if syst in histname: 
        contain_syst = True
        break
    if not contain_syst:
      if histname not in histnames_list_wosyst:
        histnames_list_wosyst.append(histname)
  histnames_list_wosyst = [h for h in histnames_list_wosyst if 'CutFlow' not in h]
  histnames_list_wosyst = [h for h in histnames_list_wosyst if 'gen' not in h]
  for sys in syl+['']:
    for histname_wosyst in histnames_list_wosyst:
      if sys != '':
        histname = histname_wosyst + '_' + sys
      else:
        histname = histname_wosyst
      for idx, samplename in enumerate(sl):
        print('fs: ', fs)
        sample = fs.Get(f'evt/tight_lep/{msoftdrop}/{samplename}')
        if not sample:
          print('sample: ', samplename)
          assert(sample)
        hist = sample.Get(histname)
        if not hist:
          histname_new = histname.replace(f'_{sys}', '')
          hist = sd.Get(histname_new)
          hist.SetName(histname)
          assert(hist)
        print('hist: ', hist.Integral(), '\t', hist.GetName())
        if idx == 0:
          addhist = hist.Clone()
          print('msoftdrop: ', msoftdrop, '\t', 'sample: ', samplename, '\t', hist.GetName(), '\t', hist.Integral(), '\t', addhist.Integral())
        else:
          addhist.Add(hist)
          print('sample: ', samplename, '\t', hist.Integral(), '\t', addhist.Integral())
      fo.cd(f'evt/tight_lep/{msoftdrop}/nonQCD')
      addhist.Write()
      fs.cd()

  for sys in syl+['']:
    for histname_wosyst in histnames_list_wosyst:
      if sys != '':
        histname = histname_wosyst + '_' + sys
      else:
        histname = histname_wosyst
      hist_data = data_obs.Get(histname)
      if not hist_data:
        histname_new = histname.replace(f'_{sys}', '')
        hist_data = data_obs.Get(histname_new)
        hist_data.SetName(histname)
        assert(hist_data)
      print('histname: ', histname)
      hist_subtract = fo.Get(f'evt/tight_lep/{msoftdrop}/nonQCD/{histname}')
      print('hist_subtract: ', hist_subtract.Integral(), '\t', hist_data.Integral())
      hist_data.Add(hist_subtract, -1)
      print('subtracted: ', hist_data.Integral())
      fo.cd(f'evt/tight_lep/{msoftdrop}/data_QCD')
      hist_data.Write()
      print('inte: ', hist_data.Integral(), '\t', msoftdrop)
