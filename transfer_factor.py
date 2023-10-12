import ROOT as r
from collections import OrderedDict
import json
'''f_deno = r.TFile('test/hadd/output/hadd_stage1.root')
f_neo = r.TFile('test/hadd/output/hadd_stage1.root')
h_deno = f_deno.Get('evt/QCD/fake_lepton/msoftdrop_L_90/HT')
h_neo = f_neo.Get('evt/QCD/fake_lepton/msoftdrop_GE_90/HT')'''
f_deno = r.TFile('/afs/cern.ch/work/s/snandan/public/myforkhtoaa/htoaa-1/test_data_ele/addBackgrounds.root')#('test.root')
f_neo = f_deno#r.TFile('test.root')
h_deno = f_deno.Get('evt/fake_lepton/msoftdrop_L_90/data_QCD/HT')
h_neo = f_neo.Get('evt/fake_lepton/msoftdrop_GE_90/data_QCD/HT')
d = OrderedDict([
    ('150', h_neo.GetBinContent(1)/h_deno.GetBinContent(1)),
    ('200', h_neo.GetBinContent(2)/h_deno.GetBinContent(2)),
    ('225', h_neo.GetBinContent(3)/h_deno.GetBinContent(3)),
    ('250', h_neo.GetBinContent(4)/h_deno.GetBinContent(4)),
    ('275', h_neo.GetBinContent(5)/h_deno.GetBinContent(5)),
    ('300', h_neo.GetBinContent(6)/h_deno.GetBinContent(6)),
    ('325', h_neo.GetBinContent(7)/h_deno.GetBinContent(7)),
    ('350', h_neo.GetBinContent(8)/h_deno.GetBinContent(8)),
    ('375', h_neo.GetBinContent(9)/h_deno.GetBinContent(9)),
    ('400', h_neo.GetBinContent(10)/h_deno.GetBinContent(10)),
    ('425', h_neo.GetBinContent(11)/h_deno.GetBinContent(11)),
    ('450', h_neo.GetBinContent(12)/h_deno.GetBinContent(12)),
    ('475', h_neo.GetBinContent(13)/h_deno.GetBinContent(13)),
    ('500', h_neo.GetBinContent(14)/h_deno.GetBinContent(14)),
    ('550', h_neo.GetBinContent(15)/h_deno.GetBinContent(15)),
    ('600', h_neo.GetBinContent(16)/h_deno.GetBinContent(16)),
    ('1000', h_neo.GetBinContent(17)/h_deno.GetBinContent(17))
    #('90', h_neo.GetBinContent(1)/h_deno.GetBinContent(1)),
    #('95', h_neo.GetBinContent(2)/h_deno.GetBinContent(2)),
    #('100', h_neo.GetBinContent(3)/h_deno.GetBinContent(3)),
    #('105', h_neo.GetBinContent(4)/h_deno.GetBinContent(4)),
    #('110', h_neo.GetBinContent(5)/h_deno.GetBinContent(5)),
    #('115', h_neo.GetBinContent(6)/h_deno.GetBinContent(6)),
    #('120', h_neo.GetBinContent(7)/h_deno.GetBinContent(7)),
    #('125', h_neo.GetBinContent(8)/h_deno.GetBinContent(8)),
    #('130', h_neo.GetBinContent(9)/h_deno.GetBinContent(9)),
    #('140', h_neo.GetBinContent(10)/h_deno.GetBinContent(10)),
    #('150', h_neo.GetBinContent(11)/h_deno.GetBinContent(11)),
    #('170', h_neo.GetBinContent(12)/h_deno.GetBinContent(12)),
    #('200', h_neo.GetBinContent(13)/h_deno.GetBinContent(13))'''
    
])

with open('transfer_factor.json', 'w') as f:
    json.dump(d, f, indent=4)
