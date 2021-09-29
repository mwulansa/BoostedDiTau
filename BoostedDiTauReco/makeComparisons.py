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
    },
    'TCP_m_50':{
        'htj_0to100':9.647e-05,
        'htj_100to400':8.817e-06,
        'htj_400toInf':1.701e-07,
    }
}

hists0={}
hists1={}
hists2={}
hists3={}
hists4={}

for sample in Samples:
    for mass in list(xsecs[sample]):
        #filename = 'h_trigStudy_Sel2_{}_w_1_{}.root'.format(sample, mass)
        filename = 'h_channelStudy_{}_w_1_{}.root'.format(sample, mass)

        fil = ROOT.TFile(filename)

        #fil2 = ROOT.TFile(filename.replace('Sel1','Sel2'))

        h0 = fil.Get('JetPt')
        
        #h1 = fil.Get('EMu_M_JetHT')
        #h2 = fil.Get('EMu_M_SingleE')
        #h3 = fil.Get('EMu_M_SingleMu')
        #h4 = fil.Get('EMu_M_MuonEG')

        #h1 = fil.Get('ETau_Cleaned_M_JetHT')
        #h2 = fil.Get('ETau_Cleaned_M_SingleE')
        #h3 = fil.Get('ETau_Cleaned_M_Tau')

        #h1 = fil.Get('TauTau_M_JetHT')
        #h2 = fil.Get('TauTau_M_Tau')

        #hname = 'ETau_Cleaned_M_Tau'
        hname = 'ETau_M_Sel2'

        h1 = fil.Get(hname)
        #h2 = fil.Get(hname.replace('Cleaned', 'Reco_Cleaned'))
        h2 = fil.Get(hname.replace('Sel2','Sel1'))

        #h1 = fil.Get('EMu_M_SingleEonly')
        #h2 = fil.Get('EMu_M_SingleMuonly')
        #h3 = fil.Get('EMu_M_MuonEGonly')

        norm=xsecs[sample][mass]*1000000*41000/fil.Get('NEvent').GetBinContent(2)
        h0.Scale(norm)
        h1.Scale(norm)

        #norm2 = xsecs[sample][mass]*1000000*41000/fil2.Get('NEvent').GetBinContent(2)
        h2.Scale(norm)
        
        #h3.Scale(norm)
        #h4.Scale(norm)
        
        if list(xsecs[sample]).index(mass)==0:
            hists0[sample] = h0
            hists1[sample] = h1
            hists2[sample] = h2
            #hists3[sample] = h3
            #hists4[sample] = h4
            hists0[sample].SetDirectory(0)
            hists1[sample].SetDirectory(0)
            hists2[sample].SetDirectory(0)
            #hists3[sample].SetDirectory(0)
            #hists4[sample].SetDirectory(0)
        else:
            hists0[sample].Add(h0)
            hists1[sample].Add(h1)
            hists2[sample].Add(h2)
            #hists3[sample].Add(h3)
            #hists4[sample].Add(h4)

Sample = 'TCP_m_30'

c1 = ROOT.TCanvas("ZLL", "ZLL", 0, 0, 800, 800)
c1.cd()
c1.SetLogy()
hists0[Sample].Draw()
c1.SaveAs(Sample+'_JetPT_trigStudy.pdf')

c2 = ROOT.TCanvas("ZLL", "ZLL", 0, 0, 800, 800)
c2.cd()
hists2[Sample].Draw('')
hists1[Sample].Draw('same')
#hists2[Sample].Draw('same')
#hists3[Sample].Draw('same')


SetHistStyle(hists1[Sample], 1, 1, 2, 0)
SetHistStyle(hists2[Sample], 6, 1, 2, 0)
#SetHistStyle(hists3[Sample], 8, 1, 2, 0)
#SetHistStyle(hists4[Sample], 9, 1, 2, 0)

l=ROOT.TLegend(0.35,0.6,0.85,0.9, "")
l.SetTextSize(0.03)
l.SetFillStyle(0)

l.AddEntry(hists1[Sample], 'Veto '+str(int(hists1[Sample].Integral())), 'alp')
l.AddEntry(hists2[Sample], 'Sequential '+str(int(hists2[Sample].Integral())), 'alp')


#l.AddEntry(hists1[Sample], 'JetHT '+str(int(hists1[Sample].Integral())), 'alp')
#l.AddEntry(hists2[Sample], 'JetHT+Tau '+str(int(hists2[Sample].Integral())), 'alp')
#l.AddEntry(hists3[Sample], 'JetHT+SingleE+Tau '+str(int(hists3[Sample].Integral())), 'alp')
#l.AddEntry(hists4[Sample], 'JetHT+SingleE+SingleMu+MuonEG '+str(int(hists4[Sample].Integral())), 'alp')

#l.AddEntry(hists1[Sample], 'EMu SingleEOnly '+str(int(hists1[Sample].Integral())), 'alp')
#l.AddEntry(hists2[Sample], 'EMu SingleMuOnly '+str(int(hists2[Sample].Integral())), 'alp')
#l.AddEntry(hists3[Sample], 'EMu MuonEGOnly '+str(int(hists3[Sample].Integral())), 'alp')

l.Draw('same')

c2.SaveAs(Sample+'_ETau_catStudy.pdf')

