from ROOT import TFile, TH1F, TH1D, TCanvas, TLegend
import os
import argparse

parser = argparse.ArgumentParser(description='htoaa analysis wrapper')
parser.add_argument('-p', dest='input_path', type=str, required=True)
parser.add_argument('-o', dest='output_file', type=str, required=True)
args=parser.parse_args()
p = args.input_path
output_file = args.output_file
inputpath = args.input_path

fori = TFile(os.path.join(inputpath, 'hadd', 'output', 'hadd_stage1.root'))
fest = fori#TFile(os.path.join(inputpath, 'met_L20', 'hadd', 'output', 'hadd_stage1_Tight_met_L20.root'))
hori = fori.Get('evt/QCD/tight_lep/msoftdrop_GE_90/hLeadingFatJetParticleNetMD_Xbb').Rebin(10)
hest = fest.Get('evt/QCD/tight_lep/msoftdrop_L_90/hLeadingFatJetParticleNetMD_Xbb')
hest.Rebin(10)
hest.SetLineColor(2)
hori.SetStats(0)
hest.SetStats(0)
print(hori, '\t', hest)
legend1 = TLegend(0.1800, 0.70, 0.90, 0.94)
legend1.SetFillStyle(0)
legend1.SetBorderSize(0)
legend1.SetFillColor(10)
legend1.SetTextSize(0.030)
legend1.AddEntry(hori, 'oroginal:{:.2f}'.format(hori.Integral()), "f")
legend1.AddEntry(hest, 'estimated:{:.2f}'.format(hest.Integral()), "f")
c = TCanvas("","",600,800)
c.SetLogy()
hori.Draw()
hest.Draw("same")
legend1.Draw('same')
c.SaveAs(output_file)
