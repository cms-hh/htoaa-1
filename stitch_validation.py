#import uproot3 as uproot
#import numpy as np
#import awkward as ak
import sys
import json
import glob
import os
import ROOT as r

f2=r.TFile('test_stitch_withcut_allsamples/hadd/output/hadd_stage1.root')#wjets_all.root')#hadd_apply_stitch.root')
f1=r.TFile('test_stitch_withcut/hadd/output/hadd_stage1.root')#wjets_inclusive/hadd/output/hadd_stage1.root')#hadd_apply_stitch_nostitch_inclusive.root')

h1=f1.Get('evt/tight_lep/msoftdrop_GE_90/WJets/gen/LHE_Njet_gen')
h2=f2.Get('evt/tight_lep/msoftdrop_GE_90/WJets/gen/LHE_Njet_gen')
legend1 = r.TLegend(0.1800, 0.75, 0.70, 0.94)

def create_canvas():
  canvas = r.TCanvas("canvas", "canvas", 800, 900);
  canvas.SetFillColor(10);
  canvas.SetFillStyle(4000)
  canvas.SetTicky(0)
  canvas.SetTickx(0)
  canvas.SetBorderSize(2);
  canvas.SetLeftMargin(0.12);
  canvas.SetBottomMargin(0.12)
  return canvas

def create_Pad(y1=0, y2=1, topmargin=0.055, leftmargin=0.12, bottommargin=0.03, rightmargin=0.12, logy=True):
  topPad = r.TPad("topPad", "topPad", 0.00, y1, 1.00, y2);
  topPad.SetFillColor(10);
  topPad.SetTopMargin(topmargin)
  topPad.SetLeftMargin(leftmargin)
  topPad.SetBottomMargin(bottommargin)
  topPad.SetRightMargin(rightmargin)
  topPad.SetLogy(logy)
  return topPad

canvas = create_canvas()
topPad =    create_Pad(y1=0.35, y2=1.0, topmargin=0.055, leftmargin=0.12, bottommargin=0.03, rightmargin=0.12, logy=True)
bottomPad = create_Pad(y1=0.0,  y2=0.35,topmargin=0.02,  leftmargin=0.12, bottommargin=0.31, rightmargin=0.12, logy=False)
canvas.cd()
topPad.Draw()
topPad.cd()


#c=r.TCanvas('','',600,800)
#c.SetLogy()
h2.SetLineColor(r.kRed)
h1.Rebin(100)
h2.Rebin(100)
h1.SetStats(0)
h2.SetStats(0)
h1.SetMarkerStyle(20)
#h1.Scale(1/h1.Integral())
#h2.Scale(1/h2.Integral())
h1.Draw('e1')
h2.Draw('hist same')
u = h1.Clone('un')
u.SetLineColor(12)
u.SetFillColorAlpha(12, 0.40)
u.Draw('e2 same')
legend1.AddEntry(h1, 'inclusive:{0:.3f}'.format(h1.Integral()), "f")
legend1.AddEntry(h2, 'inclusive+exclusive stitched:{0:.3f}'.format(h2.Integral()), "f")
legend1.AddEntry(u, 'stat Uncer', "f")
legend1.Draw('same')
canvas.cd()
bottomPad.Draw()
bottomPad.cd()
ratio = h1.Clone('ratio')
ratio.Divide(h2)
ratio.SetStats(0)
ratio.SetTitle('')
err_hist = h1.Clone("err")
err_hist.SetMarkerStyle(0)
err_hist.GetYaxis().SetTitle('inclusive / stitched')
#err_hist.GetYaxis().SetTitleSize(0.1)
for ibin in range(1, err_hist.GetNbinsX()+1):
  err_hist.SetBinContent(ibin,1)
  binerr = 0
  if h1.GetBinContent(ibin) > 0.:
    binerr = h1.GetBinError(ibin)/h1.GetBinContent(ibin)
  err_hist.SetBinError(ibin, binerr)
err_hist.GetYaxis().SetRangeUser(0.6,1.40)
err_hist.SetFillColorAlpha(12, 0.40)
err_hist.Draw('same,e2')
line = r.TF1("line","1", ratio.GetBinLowEdge(1), ratio.GetBinLowEdge(ratio.GetNbinsX())+ratio.GetBinWidth(ratio.GetNbinsX()));
line.SetLineStyle(3);
line.SetLineColor(r.kBlack);
line.Draw("same");
ratio.Draw('e1psame')
canvas.SaveAs('stitch_HT.png')
