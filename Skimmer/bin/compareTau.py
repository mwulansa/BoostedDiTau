import os, sys, math
import ROOT
#from myPlotFunctions import *
import numpy as np

mass="50"

fname1 = "/uscms/home/jingyu/nobackup/TCP/boostedDiTauReco/CMSSW_8_0_30/src/BoostedDiTau/BoostedDiTauReco/h_plotBoostedDiTauReco_m"+mass+"_v4_80X_backup.root"

fname2 = "h_plot_tauIds_m"+mass+"_Medium.root"

ROOT.gROOT.ForceStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetPadTopMargin(0.05)
ROOT.gStyle.SetPadRightMargin(0.04)
ROOT.gStyle.SetPadLeftMargin(0.09)

#histName = "hMuTau_MuCleaned_MuPt"
#histName = "hNMuons"
histName1 = "h_tauPt"
histName2 = "h_tauPt_new"


#histName = "hETau_ECleaned_EPt"
#histName = "hETau_ECleaned_EPt_genMatched"
#histName = "hMuTau_M"
#histName = "hETau_ECleaned_EPt"
#histName = "h_isoEPt"

f1=ROOT.TFile(fname1)
f2=ROOT.TFile(fname2)

h1 = f1.Get(histName1)

h2 = f2.Get(histName1)

h3 = f2.Get(histName2)

h1.Rebin(5)
h2.Rebin(5)
h3.Rebin(5)


c=ROOT.TCanvas("myCanvas", "myCanvas")

c.SetLogy()

h1.SetLineColor(2)

h2.Draw()
#h1.SetMinimum(2)
h1.Draw("same")
h2.SetLineColor(ROOT.kGreen+3)
h3.SetLineColor(4)
h3.Draw("same")

l=ROOT.TLegend(0.7, 0.6, 0.9, 0.9, "")
l.SetBorderSize(0)
l.AddEntry(h1, "80X", "alp")
l.AddEntry(h2, "94X oldID", "alp")
l.AddEntry(h3, "94X newID", "alp")
l.SetTextSize(0.04)


l.Draw("same")

c.SaveAs("comp_tauId_m"+mass+"_Medium.pdf")
