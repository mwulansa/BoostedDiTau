import os, sys, math
import ROOT
from myPlotFunctions import *
import numpy as np

mass = "50"

fname1 = "/uscms/home/jingyu/nobackup/TCP/boostedDiTauReco/CMSSW_8_0_30/src/BoostedDiTau/BoostedDiTauReco/h_plotBoostedDiTauReco_m"+mass+"_v4_80X_backup.root"

#fname1 = "/uscms/home/jingyu/nobackup/TCP/boostedDiTauReco/CMSSW_8_0_30/src/BoostedDiTau/BoostedDiTauReco/h_plotBoostedDiTauReco_m10_v4_80X_backup.root"

fname2 = "h_plotBoostedDiTauReco_m"+mass+"_v4_94X_backup.root"

#fname2 = "h_plotBoostedDiTauReco_m10_v4_94X_fullStats.root"

#fname2 = "h_plotBoostedDiTauReco_m10_v4_94X_noSC.root"

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
#histName = "h_tauPt"
#histName = "hETau_ECleaned_M"
#histName = "hETau_ECleaned_EPt"
#histName = "hETau_ECleaned_EPt_genMatched"
#histName = "hMuTau_M"
histName = "hETau_ECleaned_EPt"
#histName = "hMuTau_MuCleaned_dR"
#histName = "h_isoEPt"

f1=ROOT.TFile(fname1)
f2=ROOT.TFile(fname2)

h1 = f1.Get(histName)

h2 = f2.Get(histName)

h1.Rebin(5)
h2.Rebin(5)

#p1,p2,p3,p4,p5,p6,p7=RatioPlotTopXaxisSetup
#SetupAxis(h1, p1,p2,p3,p4,p5,p6,p7)
#h4.GetYaxis().SetMaxDigits(2)
#p1,p2,p3,p4,p5,p6,p7=RatioPlotTopYaxisSetup
#p5='N'
#SetupAxis(h1, p1,p2,p3,p4,p5,p6,p7)

c = ROOT.TCanvas("Comparisons", "Comparisons", 0, 0, 650, 800)
c.Divide(1,1)
c.cd()
pad1=ROOT.TPad("","",0,0,1,0.5)
pad1.SetTopMargin(0.08)
pad1.SetBottomMargin(0.25)
c.cd()
pad2=ROOT.TPad("","",0,0.5,1,1)
pad2.SetBottomMargin(0.001)
pad2.SetTopMargin(0.07)

pad2.SetLogy()

pad2.Draw()
pad2.cd()

h1.SetLineColor(2)
h1.Draw()
#h1.SetMinimum(2)
h2.Draw("same")
h2.SetLineColor(ROOT.kGreen+3)

l=ROOT.TLegend(0.6, 0.7, 0.8, 0.9, "")
SetupLegend(l, 0.04)
l.AddEntry(h1, "80X", "alp")
l.AddEntry(h2, "94X", "alp")
l.SetBorderSize(0)
l.Draw("same")

c.cd()
pad1.Draw()
pad1.cd()
h2_clone = h2.Clone("den")
h1_clone = h1.Clone("num")
h2_clone.Divide(h1_clone)
h1_clone.Divide(h1_clone)

h2_clone.GetYaxis().CenterTitle()
h2_clone.GetYaxis().SetTitle("")
h2_clone.SetMaximum(2)
h2_clone.SetMinimum(0)
#SetRatioPlotTick(h2_clone)

h2_clone.Draw()
h1_clone.Draw("same")

c.SaveAs('comp_m'+mass+"_"+histName+'.pdf')
