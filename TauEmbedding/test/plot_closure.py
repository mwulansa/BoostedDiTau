import os, sys, math
import ROOT

from myPlotFunctions import *

ROOT.gROOT.ForceStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetPadLeftMargin(0.13)
ROOT.gStyle.SetPadBottomMargin(0.11)
ROOT.gStyle.SetPadTopMargin(0.09)
ROOT.gStyle.SetPadRightMargin(0.07)
ROOT.gStyle.SetLegendBorderSize(0)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)

ROOT.gStyle.SetOptFit(1)
#ROOT.gStyle.SetPaintTextFormat("4.1f")

whichVar1 = 'mtautau_deep_cor'
whichVar2 = 'mtautauEmbed_deep_cor'
#whichVar1 = 'tau1Pt_deep'
#whichVar2 = 'tau1PtEmbed_deep_cor'
#whichVar1 = 'genTau1Pt'
#whichVar2 = 'genTau1Pt'

#whichVar1 = 'genTauDR'
#whichVar2 = 'genTauDR'

#whichVar1 = 'dRtautau_deep_cor'
#whichVar2 = 'dRtautauEmbed_deep_cor'
#whichVar1 = 'mu1Pt_cor'
#whichVar2 = 'mu1PtEmbed_cor'
#whichVar1 = 'ePt'
#whichVar2 = 'ePtEmbed'

#whichVar1 = 'mmumu_med'
#whichVar2 = 'mmumuEmbed_med'

#whichVar1 = 'dRmumu_cor'
#whichVar2 = 'dRmumuEmbed_cor'

#whichVar1 = 'memu_cor'
#whichVar2 = 'memuEmbed_cor'
#whichVar1 = 'muPt'
#whichVar2 = 'tauPtEmbed'

#whichVar1 = 'mmumu'
#whichVar2 = 'mtautauEmbed_deep'

#whichVar1 = 'mu2Pt_cor'
#whichVar2 = 'mu2PtEmbed_cor'

#whichVar1 = 'genTaumu2Pt'
#whichVar2 = 'genTaumu2Pt'

#whichVar1 = 'genTauM'
#whichVar2 = 'mmumuEmbed_genM'

rebin = 1

f = ROOT.TFile('h_boosted.root')
#f = ROOT.TFile('h_muon_embedded_boosted.root')
#f = ROOT.TFile('h_embedded_boosted.root')
h = f.Get(whichVar1)
#h.Add(f.Get(whichVar2))
h.Scale(500/497)
h.Rebin(rebin)

print(h.Integral())

fembed = ROOT.TFile('h_embedded_boosted.root')
#fembed = ROOT.TFile('h_muon_embedded_boosted.root')
hembed = fembed.Get(whichVar2)
#hembed.Add(fembed.Get(whichVar2))
#hembed = f.Get(whichVar2)
#hembed.Scale(500/496)
hembed.Scale(500000/435889)
hembed.Rebin(rebin)

print(hembed.Integral())

l=ROOT.TLegend(0.55,0.6,0.85,0.9, "")
l.SetTextSize(0.04)
l.SetFillStyle(0)
#SetupLegend(l, 0.04)



c = ROOT.TCanvas("ZLL", "ZLL", 0, 0, 650, 800)
c.Divide(1,1)
c.cd()
padDown=ROOT.TPad("","",0,0,1,0.3)
padDown.SetTopMargin(0.08)
padDown.SetBottomMargin(0.25)
#padDown.Draw()
c.cd()
padUp=ROOT.TPad("","",0,0.3,1,1)
padUp.SetBottomMargin(0.001)
padUp.SetTopMargin(0.07)
#padUp.Draw()

#c.SetLogy()

padUp.Draw()
padUp.cd()
#padUp.SetLogy()
#h.Draw('hist')
h.SetLineColor(2)
#h.SetMinimum(100)
hembed.Draw('hist')
hembed.SetMarkerSize(1)
hembed.SetMarkerStyle(20)
hembed.SetLineColor(1)
h.Draw('hist')
hembed.Draw('samehist')
l.Draw('same')
#bg.SetMinimum(1000)
p1,p2,p3,p4,p5,p6,p7=RatioPlotTopXaxisSetup
SetupAxis(h, p1,p2,p3,p4,p5,p6,p7)
#h4.GetYaxis().SetMaxDigits(2)
p1,p2,p3,p4,p5,p6,p7=RatioPlotTopYaxisSetup
p5='N'
SetupAxis(h, p1,p2,p3,p4,p5,p6,p7)

c.cd()
padDown.Draw()
padDown.cd()
#padDown.cd()
ratio = hembed.Clone("Ratio")
#h_bg = bg.Sum()
ratio.Divide(h)
#ratio = h.Clone("Ratio")
#ratio.Divide(hembed)



function=ROOT.TF1("Fit_correction","pol1",50,500)
function.SetLineColor(ROOT.kGreen)
function.SetParameters(0,1.)
function.SetParameters(0,0.)
ratio.Fit(function,"0R")
ratio.Draw('')
#ratio.SetMarkerSize(1.8)
function.Draw("same")

#ratio.GetXaxis().SetTickSize(0.4)
ratio.SetMaximum(2)
ratio.SetMinimum(0.)
p1,p2,p3,p4,p5,p6,p7=RatioPlotBottomXaxisSetup
## if whichVar == 'mll':
##     p5='m_{ll}'
## elif whichVar == 'met':
##     p5='MET'
## elif whichVar == 'dRll':
##     p5 = 'dR_{ll}'
## elif whichVar == 'jetPt':
##     p5 = 'p_{t,jet}'
SetupAxis(ratio, p1,p2,p3,p4,p5,p6,p7)
p1,p2,p3,p4,p5,p6,p7=RatioPlotBottomYaxisSetup
p5='Ratio'
SetupAxis(ratio, p1,p2,p3,p4,p5,p6,p7)
ratio.GetYaxis().CenterTitle()

#data.Draw()
#data.GetXaxis().SetRangeUser(0,100)
#data.SetMinimum(10000)
#hists['DYJetsToLL'].Draw('h')
#hists['DYJetsToLL'].GetXaxis().SetRangeUser(0,200)
#hists['DYJetsToLL'].SetFillColor(ROOT.kOrange-3)
#hists['DYJetsToLL'].SetLineColor(ROOT.kOrange-3)
#data.Draw('same')
#c.SaveAs('closure_muonEmbed_{}.pdf'.format(whichVar1))
c.SaveAs('closure_{}.pdf'.format(whichVar1))
