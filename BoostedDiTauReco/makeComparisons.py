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

Samples=['TCP_m_30']

xsecs={
    'TCP_m_30':{
        'htj_0to100':9.674e-05,
        'htj_100to400':9.457e-06,
        'htj_400toInf':1.623e-07,
    }
}

hists0={}
hists1={}
hists2={}
hists3={}
hists4={}

for sample in Samples:
    for mass in list(xsecs[sample]):
        filename = 'h_trigStudy_{}_w_1_{}.root'.format(sample, mass)

        fil = ROOT.TFile(filename)

        h0 = fil.Get('JetPt')
        h1 = fil.Get('EMu_M_JetHT')
        h2 = fil.Get('EMu_M_SingleE')
        h3 = fil.Get('EMu_M_SingleMu')
        h4 = fil.Get('EMu_M_MuonEG')

        #h1 = fil.Get('ETau_M_JetHT')
        #h2 = fil.Get('ETau_M_SingleE')
        #h3 = fil.Get('ETau_M_Tau')

        #h1 = fil.Get('TauTau_M_JetHT')
        #h2 = fil.Get('TauTau_M_Tau')

        norm=xsecs[sample][mass]*1000000*41000/fil.Get('NEvent').GetBinContent(2)
        h0.Scale(norm)
        h1.Scale(norm)
        h2.Scale(norm)
        h3.Scale(norm)
        h4.Scale(norm)
        
        if list(xsecs[sample]).index(mass)==0:
            hists0[sample] = h0
            hists1[sample] = h1
            hists2[sample] = h2
            hists3[sample] = h3
            hists4[sample] = h4
            hists0[sample].SetDirectory(0)
            hists1[sample].SetDirectory(0)
            hists2[sample].SetDirectory(0)
            hists3[sample].SetDirectory(0)
            hists4[sample].SetDirectory(0)
        else:
            hists0[sample].Add(h0)
            hists1[sample].Add(h1)
            hists2[sample].Add(h2)
            hists3[sample].Add(h3)
            hists4[sample].Add(h4)

c1 = ROOT.TCanvas("ZLL", "ZLL", 0, 0, 800, 800)
c1.cd()
c1.SetLogy()
hists0['TCP_m_30'].Draw()
c1.SaveAs(sample+'_JetPT_trigStudy.pdf')

c2 = ROOT.TCanvas("ZLL", "ZLL", 0, 0, 800, 800)
c2.cd()
hists4['TCP_m_30'].Draw('')
hists1['TCP_m_30'].Draw('same')
hists2['TCP_m_30'].Draw('same')
hists3['TCP_m_30'].Draw('same')


SetHistStyle(hists1['TCP_m_30'], 1, 1, 2, 0)
SetHistStyle(hists2['TCP_m_30'], 6, 1, 2, 0)
SetHistStyle(hists3['TCP_m_30'], 8, 1, 2, 0)
SetHistStyle(hists4['TCP_m_30'], 9, 1, 2, 0)

l=ROOT.TLegend(0.35,0.6,0.85,0.9, "")
l.SetTextSize(0.03)
l.SetFillStyle(0)

l.AddEntry(hists1['TCP_m_30'], 'JetHT '+str(int(hists1['TCP_m_30'].Integral())), 'alp')
l.AddEntry(hists2['TCP_m_30'], 'JetHT+SingleE '+str(int(hists2['TCP_m_30'].Integral())), 'alp')
l.AddEntry(hists3['TCP_m_30'], 'JetHT+SingleE+SingleMu '+str(int(hists3['TCP_m_30'].Integral())), 'alp')
l.AddEntry(hists4['TCP_m_30'], 'JetHT+SingleE+SingleMu+MuonEG '+str(int(hists4['TCP_m_30'].Integral())), 'alp')

l.Draw('same')

c2.SaveAs(sample+'_EMu_trigStudy.pdf')

