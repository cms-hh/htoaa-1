import ROOT as r
import collections
import sys
dprocs = collections.OrderedDict()
dprocs["data_fakes"]   = {"color" :  1, "fillStype" : 3345, "label" : "Fakes"  ,    "make border" : True}
dprocs["Fakes"]   = {"color" :  1, "fillStype" : 3345, "label" : "Fakes"  ,    "make border" : True}
dprocs["Convs"]      = {"color" :   800, "fillStype" : 1001, "label" : "Conversions", "make border" :  True}
dprocs["TT"]         = {"color" : 822, "fillStype" : 1001, "label" : 't#bar{t} + jets'   , "make border" : True
}
dprocs["ST"]         = {"color" : 823, "fillStype" : 1001, "label" : 'ST'   , "make border" : True}
dprocs["data_QCD"]      = {"color" : 851, "fillStype" : 1001, "label" : "QCD"       , "make border":True}
dprocs["Other_bbWW"]      = {"color" : 851, "fillStype" : 1001, "label" : "Other_bbWW"       , "make border":True}
dprocs["WZ"]          = {"color" : 7, "fillStype" : 1001, "label" : "WZ"         , "make border" :True}
dprocs["WJets"]          = {"color" : 634, "fillStype" : 1001, "label" : "WJets"         , "make border" :True}
dprocs["WZZ"]         = {"color" : 16,   "fillStype" : 1001, "label" : "WZZ"          , "make border" : True}
dprocs["ZZZ"]        = {"color" : 822, "fillStype" : 1001, "label" : "ZZZ"          , "make border" : True}
dprocs["ZZ"]        = {"color" : 5, "fillStype" : 1001, "label" : "ZZ"        , "make border" : False}
dprocs["ttH"]       = {"color" : 4, "fillStype" : 1001, "label" : "t#bar{t}H"  , "make border" : True}
dprocs["Z"]         = {"color" : 610, "fillStype" : 1001, "label" : "ZtoQQ"         , "make border" : True}
dprocs["SH"]         = {"color" : 628, "fillStype" : 1001, "label" : "Single H"         , "make border" : False
}

def set_axis(axis, title='', titleoffset=1.25, labelsize=10, labeloffset=1000, labelcolor=10, titlecolor=10):
   axis.SetTitle(title)
   axis.SetTitleOffset(titleoffset)
   axis.SetLabelSize(labelsize)
   axis.SetLabelOffset(labeloffset)
   axis.SetLabelColor(labelcolor)
   axis.SetTitleColor(titlecolor)

def create_legend():
  legend_y0 = 0.75
  legendPosX = 0.570
  legendPosY = 0.510
  legendSizeX = 0.360
  legendSizeY = 0.420
  legend1 = r.TLegend(0.1800, legend_y0, 0.90, 0.94)
  legend1.SetNColumns(3)
  legend1.SetFillStyle(0)
  legend1.SetBorderSize(0)
  legend1.SetFillColor(10)
  legend1.SetTextSize(0.030)
  return legend1

def set_histproperty(hist_rebin_local, itemDict, addlegend=True):
  #hist_rebin_local.SetMarkerSize(0)
  hist_rebin_local.SetFillColor(itemDict["color"])
  if not itemDict["fillStype"] == 0 :
    hist_rebin_local.SetFillStyle(itemDict["fillStype"])
  if "none" not in itemDict["label"] and addlegend :
        legend1.AddEntry(hist_rebin_local, itemDict["label"], "f")
  if itemDict["make border"] == True :
    hist_rebin_local.SetLineColor(1)
  else :
    hist_rebin_local.SetLineColor(itemDict["color"])
  hist_rebin_local.GetXaxis().SetTickLength(0.04)
  hist_rebin_local.GetYaxis().SetTitleSize(0.055)
  hist_rebin_local.GetYaxis().SetLabelSize(0.050)

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

f = r.TFile('test_data_ele_tight_lepton/hadd.root')
d = f.Get('evt/tight_lep/msoftdrop_GE_90')
#d = uproot.open('test/hadd/output/hadd_stage1.root')['evt']
bkg_hists = []
legend1 = create_legend()
histogramStack_mc = r.THStack()
var = 'hLeadingFatJetParticleNetMD_Xbb'
for idx, key in enumerate([k.GetName() for k in d.GetListOfKeys() if 'nonQCD' not in k.GetName()]):
   dir = d.Get(key)
   print(dir)
   #hist = dir.Get('tight_lep/msoftdrop_GE_90/'+var) if 'data_QCD' not in key\
   #else dir.Get('tight_lep/msoftdrop_L_90/'+var)
   hist = dir.Get(var) if 'data_QCD' not in key\
          else f.Get('evt/tight_lep/msoftdrop_L_90/data_QCD/'+var)
   hist.SetName(key)
   if 'data_obs' in key:
      data_hist = hist#dir.Get(var)
      print(data_hist.Integral())
      data_hist.Rebin(10)
      continue
   if 'SUSY' in key: continue
   hist.Rebin(10)
   if idx == 0:
      total_hist = hist.Clone('total_hist')
   else:
      total_hist.Add(hist)
   hist.SetMarkerSize(0)
   set_histproperty(hist, dprocs[key], True)
   bkg_hists.append(hist)

bkg_hists = sorted(bkg_hists, key=lambda x: x.Integral())
for hist in bkg_hists:
   print(hist.GetName(), '\t', hist.Integral())
   histogramStack_mc.Add(hist)

canvas = create_canvas()
topPad =    create_Pad(y1=0.35, y2=1.0, topmargin=0.055, leftmargin=0.12, bottommargin=0.03, rightmargin=0.12, logy=True)
bottomPad = create_Pad(y1=0.0,  y2=0.35,topmargin=0.02,  leftmargin=0.12, bottommargin=0.31, rightmargin=0.12, logy=False)
canvas.cd()
topPad.Draw()
topPad.cd()
xAxis_top = total_hist.GetXaxis();
set_axis(xAxis_top, title='LeadingFatJetParticleNetMD_Xbb', titleoffset=1.25, labelsize=10, labeloffset=1000, labelcolor=10, titlecolor=10)
yAxis_top = total_hist.GetYaxis()
#set_axis(xAxis_top, title='DNN Score in %s' %evtcat, titleoffset=0.70, labelsize=0.05, labeloffset=1000)                                    
yAxis_top.SetTitle("dN/d(var)")
yAxis_top.SetTitleOffset(0.70)
yAxis_top.SetTitleSize(0.06)
yAxis_top.SetLabelSize(0.05)
yAxis_top.SetTickLength(0.055)

max_ = 100*total_hist.GetMaximum()
total_hist.GetYaxis().SetRangeUser(0.1, max_)
total_hist.SetStats(False)
total_hist.SetTitle('')
total_hist.SetMarkerColor(16)#0)                                                                                                             
#total_hist.SetMarkerStyle(20)
total_hist.SetMarkerSize(0)
total_hist.SetLineWidth(0)
total_hist.SetFillColorAlpha(12, 0.40)
total_hist.Draw("e2,same")
total_hist.SetMarkerColor(16)
total_hist.SetMarkerSize(0)
total_hist.SetLineWidth(0)
total_hist.SetFillColorAlpha(12, 0.40)
total_hist.Draw("e2,same")
legend1.AddEntry(total_hist, "Uncertainty", "f")


histogramStack_mc.Draw("hist same")
total_hist.Draw("axis,same")
legend1.Draw("same")

data_hist.SetMarkerSize(2);
data_hist.SetMarkerColor(r.kBlack);
data_hist.SetLineColor(r.kBlack);
data_hist.SetMarkerStyle(20);
data_hist.Draw("e1P same")

canvas.cd()
bottomPad.Draw()
bottomPad.cd()
bottomPad.SetLogy(0)

err_hist = total_hist.Clone("err")

for ibin in range(1, err_hist.GetNbinsX()+1):
   err_hist.SetBinContent(ibin,0)
   if total_hist.GetBinContent(ibin) > 0.:
      binerr = total_hist.GetBinError(ibin)/total_hist.GetBinContent(ibin)
      err_hist.SetBinError(ibin, binerr)

err_hist.SetMinimum(-0.3)
err_hist.SetMaximum(0.3)
xAxis_bottom = err_hist.GetXaxis();
xAxis_bottom.SetTitle(xAxis_top.GetTitle());
xAxis_bottom.SetLabelColor(1);
xAxis_bottom.SetTitleColor(1);
xAxis_bottom.SetTitleOffset(1.20);
xAxis_bottom.SetTitleSize(0.12);
xAxis_bottom.SetLabelOffset(0.02);
xAxis_bottom.SetLabelSize(0.10);
xAxis_bottom.SetTickLength(0.055);

yAxis_bottom = err_hist.GetYaxis();
yAxis_bottom.SetTitle("#frac{Data - Simulation}{Simulation}");
yAxis_bottom.SetTitleOffset(0.60);
yAxis_bottom.SetLabelSize(0.06)
yAxis_bottom.SetNdivisions(505);
yAxis_bottom.CenterTitle();
yAxis_bottom.SetTitleSize(0.09);

err_hist.Draw("e2, same")
line = r.TF1("line","0", xAxis_bottom.GetXmin(), xAxis_bottom.GetXmax());
line.SetLineStyle(3);
line.SetLineColor(r.kBlack);
line.Draw("same");

histogramRatio = data_hist.Clone("histogramRatio");
for ibin in range(0, total_hist.GetNbinsX()):
   bincontent = total_hist.GetBinContent(ibin+1)
   if bincontent == 0: continue
   if bincontent >0:
      histogramRatio.SetBinContent(ibin+1, data_hist.GetBinContent(ibin+1)/bincontent -1)
      histogramRatio.SetBinError(ibin+1, data_hist.GetBinError(ibin+1)/data_hist.GetBinContent(ibin+1))
   else:
      histogramRatio.SetBinContent(ibin+1, -1)
histogramRatio.SetMarkerStyle(data_hist.GetMarkerStyle());
histogramRatio.SetMarkerSize(data_hist.GetMarkerSize());
histogramRatio.SetMarkerColor(data_hist.GetMarkerColor());
histogramRatio.SetLineColor(data_hist.GetLineColor());
histogramRatio.Draw("e1psame")
for e in ['pdf', 'png']:
   plot = 'plot.%s' %(e)
   canvas.SaveAs(plot)
del canvas
del histogramStack_mc

'''
canvas.cd()
del topPad
del bottomPad
canvas.SaveAs('plot.png')
'''
