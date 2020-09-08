import os, sys, math
import ROOT
from myPlotFunctions import *
import numpy as np

#mass = "50"

fname1 = "/uscms/home/jingyu/nobackup/TCP/boostedDiTauReco/CMSSW_8_0_30/src/BoostedDiTau/BoostedDiTauReco/h_plot_electronIds_m_10_80X.root"

#fname1 = "/uscms/home/jingyu/nobackup/TCP/boostedDiTauReco/CMSSW_8_0_30/src/BoostedDiTau/BoostedDiTauReco/h_plotBoostedDiTauReco_m10_v4_80X_backup.root"


fname2 = "h_plot_electronIds_m_10_94X.root"
fname3 = "h_plot_electronIds_m_50_120_94X.root"

saveName = "comp_electronIds_Zee_94X.pdf"

#fname2 = "h_plotBoostedDiTauReco_m10_v4_94X_fullStats.root"

ROOT.gROOT.ForceStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetPadTopMargin(0.05)
ROOT.gStyle.SetPadRightMargin(0.04)
ROOT.gStyle.SetPadLeftMargin(0.09)

fname = fname3

histName1 = "gen_ePt"
histName2 = "ePt"
histName3 = "gen_matched_ePt"

f=ROOT.TFile(fname)

h1 = f.Get(histName1)
h2 = f.Get(histName2)
h3 = f.Get(histName3)

rebin = 10

h1.Rebin(rebin)
h2.Rebin(rebin)
h3.Rebin(rebin)

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
pad1.SetTopMargin(0.0)
pad1.SetBottomMargin(0.15)
c.cd()
pad2=ROOT.TPad("","",0,0.5,1,1)
pad2.SetBottomMargin(0.00)
pad2.SetTopMargin(0.07)

pad2.SetLogy()

pad2.Draw()
pad2.cd()

h1.SetLineColor(2)
h1.SetLineWidth(2)
h1.Draw()
#h1.SetMinimum(2)
h2.Draw("same")
h2.Sumw2()
h2.SetLineColor(ROOT.kGreen+3)
h2.SetLineWidth(2)

h3.Draw("same")
h3.SetLineColor(4)
h3.SetLineWidth(2)
h3.Sumw2()

l=ROOT.TLegend(0.6, 0.6, 0.8, 0.85, "")
SetupLegend(l, 0.04)
l.AddEntry(h1, "Gen", "alp")
l.AddEntry(h2, "Reco", "alp")
l.AddEntry(h2, "Reco Gen Matched", "alp")
l.SetBorderSize(0)
l.Draw("same")

c.cd()
pad1.Draw()
pad1.cd()
h1_clone = h1.Clone("den")
h2_clone = h2.Clone("num1")
h3_clone = h3.Clone("num2")
h2_clone.Divide(h1_clone)
h3_clone.Divide(h1_clone)
h1_clone.Divide(h1_clone)

h1_clone.GetYaxis().CenterTitle()
h1_clone.GetYaxis().SetTitle("")
h1_clone.SetMaximum(1.49)
h1_clone.SetMinimum(0)
#SetRatioPlotTick(h2_clone)

h1_clone.Draw("hist")
h2_clone.Draw("same")

h3_clone.Draw("same")

c.SaveAs(saveName)
