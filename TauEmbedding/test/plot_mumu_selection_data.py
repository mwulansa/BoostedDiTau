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

xsecs={
    'DYJetsToLL':{
        'M-1to5_HT-70to100':900.0*1.23,
        'M-1to5_HT-100to200':671.3*1.23, 
        'M-1to5_HT-200to400':118.3*1.23, 
        'M-1to5_HT-400to600':11.72*1.23, 
        'M-1to5_HT-600toInf':3.428*1.23, 
        'M-5to50_HT-70to100':303.9*1.23,
        'M-5to50_HT-100to200':220.8*1.23, 
        'M-5to50_HT-200to400':37.87*1.23, 
        'M-5to50_HT-400to600':3.628*1.23, 
        'M-5to50_HT-600toInf':1.107*1.23, 
        'M-50':4963.0*1.23
    },
    'DYJetsToLLNLO':{
        'M-10to50':18610,
        'M-50':5765.4
    },
    'TTJets':{
        'DiLept':87.315
    },
    'TTJetsNLO':{
        'DiLept':76.75
    },
    'TTTo2L2Nu':{
        'inclusive':87.31
    },
    'WJetsToLNu':{
        'inclusive':61526.7
    },
    'WJetsToLNuNLO':{
        'inclusive':61526.7
    },
    'Diboson':{
        'WW':118.7,
        'WZ':47.13,
        'ZZ':16.523
    },
    'QCDMuEnriched':{
        'Pt-15to20':3.82e+06,
        'Pt-20to30':2.96e+06,
        'Pt-30to50':1.65e+06,
        'Pt-50to80':4.38e+05,
        'Pt-80to120':1.06e+05,
        'Pt-120to170':2.52e+04,
        'Pt-170to300':8.65e+03,
        'Pt-300to470':7.97e+02,
        'Pt-470to600':7.90e+01,
        'Pt-600to800':2.51e+01,
        'Pt-800to1000':4.71,
        'Pt-1000toInf':1.62
    },
    'QCD':{
        'Pt-15to30':1.84e+9,
        'Pt-30to50':1.41e+8,
        'Pt-50to80':1.92e+7,
        'Pt-80to120':2.76e+6,
        'Pt-120to170':4.71e+5,
        'Pt-170to300':1.17e+5,
        'Pt-300to470':7.82e+3,
        'Pt-470to600':6.48e+2,
        'Pt-600to800':1.87e+2,
        'Pt-800to1000':3.23e+1,
        'Pt-1000to1400':9.41,
        'Pt-1400to1800':0.84,
        'Pt-1800to2400':0.11,
        'Pt-2400to3200':6.82e-3,
        'Pt-3200toInf':1.65e-4
    }
}


whichTree = 'mumuInfoBoosted'

whichFolder = 'analyzeMuonsForEmbeddingLooseIso'

whichVar = 'mll'
#whichVar = 'dRll'
#whichVar = 'met'
#whichVar = 'jetPt'

if whichFolder == 'analyzeMuonsForEmbeddingLooseIso':
    whichFolderGen = 'analyzeMuMuGenLooseIso'
elif whichFolder == 'analyzeMuonsForEmbeddingLooseIsoDB':
    whichFolderGen = 'analyzeMuMuGenLooseIsoDB'
else:
    whichFolderGen = 'analyzeMuMuGen'

if whichVar == 'mll':
    var = '.mass'
    binning = '(10000, 0, 100)'
elif whichVar == 'dRll':
    var = '.dR'
    binning = '(80,0,8)'
elif whichVar == 'met':
    whichTree = 'globalInfo'
    var = '.met'
    binning = '(1000, 0, 1000)'
elif whichVar == 'jetPt':
    whichTree = 'jetInfo'
    var = '.pt'
    binning = '(1000, 0, 1000)'


Samples=['WJetsToLNuNLO','QCD','Diboson','TTTo2L2Nu','DYJetsToLLNLO']

hists={}
for sample in Samples:
    for mass in list(xsecs[sample]):
        if sample == 'WJetsToLNuNLO':
            filename = 'ntuple_{}_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v3.root'.format(sample)
        elif sample == 'Diboson':
            filename = 'ntuple_{}_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3.root'.format(mass)
        elif sample == 'TTJetsNLO':
            filename = 'ntuple_{}_{}_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8_v3.root'.format(sample, mass)
        elif sample == 'DYJetsToLLNLO':
            filename = 'ntuple_{}_{}_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v3.root'.format(sample, mass)
        elif sample == 'TTTo2L2Nu':
            filename = 'ntuple_{}_TuneCP5_PSweights_13TeV-powheg-pythia8_v3.root'.format(sample)
        else:
            filename = 'ntuple_{}_{}_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3.root'.format(sample, mass)
        print(filename)
        fil = ROOT.TFile(filename)
        t_zll = fil.Get(whichFolder + '/analysisTree')
        t_gen = fil.Get(whichFolderGen+'/genTree')
        t_zll.AddFriend(t_gen, 'genInfo')
        norm = xsecs[sample][mass]*35900/fil.Get('analyzeLumiMC/Nevts').GetBinContent(2)
        print(norm)
        #z_ll.Scale(norm)

        
        t_zll.Draw(whichTree+var+' >> tt'+binning, 'puWeight*genWeight*(globalInfo.Nb < 0.5)')
        tt = ROOT.gDirectory.Get("tt")
        tt.Scale(norm)


        if 'DYJetsToLLNLO' in sample:
            t_zll.Draw(whichTree+var+' >> mumu'+binning, 'puWeight*genWeight*(genLmInfo.pdgID < 14 && globalInfo.Nb < 0.5)')
            t_zll.Draw(whichTree+var+' >> tautau'+binning, 'puWeight*genWeight*(genLmInfo.pdgID > 14 && globalInfo.Nb < 0.5)')
            mumu = ROOT.gDirectory.Get("mumu")
            tautau = ROOT.gDirectory.Get("tautau")
            mumu.Scale(norm)
            tautau.Scale(norm)
        
        if list(xsecs[sample]).index(mass)==0:
            hists[sample] = tt    
            hists[sample].SetDirectory(0)

            if 'DYJets' in sample:
                h1 = mumu.Clone(sample.replace("LL", "MM"))
                hists["DYJetsToMM"] = h1
                hists["DYJetsToMM"].SetDirectory(0)
                h2 = tautau.Clone(sample.replace("LL", "TT"))
                hists["DYJetsToTT"] = h2
                hists["DYJetsToTT"].SetDirectory(0)
        else:
            hists[sample].Add(tt)
            if 'DYJets' in sample:
                hists["DYJetsToMM"].Add(mumu)
                hists["DYJetsToTT"].Add(tautau)


fdata = ROOT.TFile("ntuple_2016_DoubleMu_v3.root")
t_data = fdata.Get(whichFolder + '/analysisTree')
t_data.Draw(whichTree+var+' >> data'+binning,'globalInfo.Nb < 0.5')
data = ROOT.gDirectory.Get("data")
#print data


bg=ROOT.THStack("bg","Stacked 1D histograms")

l=ROOT.TLegend(0.55,0.6,0.85,0.9, "")
l.SetTextSize(0.04)
l.SetFillStyle(0)
#SetupLegend(l, 0.04)

print(hists)

Samples=['WJetsToLNuNLO','Diboson','TTTo2L2Nu','DYJetsToTT','QCD','DYJetsToMM']
#Samples=['WJetsToLNuNLO','Diboson','TTJetsNLO','DYJetsToTT','DYJetsToMM']
#Samples=['WJetsToLNu','Diboson','TTJetsNLO','DYJetsToLLNLO']
rebin = 5

for Sample in Samples:
    print(Sample, hists[Sample].Integral())
    hists[Sample].Rebin(rebin)
    #hists[Sample].GetXaxis().SetRangeUser(0,100)
    
    if Sample=='DYJetsToMM':
        hists[Sample].SetFillColor(ROOT.kOrange-3)
        hists[Sample].SetLineColor(1)
    elif Sample=='DYJetsToTT':
        hists[Sample].SetFillColor(ROOT.kMagenta-2)
        hists[Sample].SetLineColor(1)
    elif Sample=='TTTo2L2Nu':
        hists[Sample].SetFillColor(ROOT.kAzure+2)
        hists[Sample].SetLineColor(1)
    elif Sample=='WJetsToLNuNLO':
        hists[Sample].SetFillColor(ROOT.kRed+2)
        hists[Sample].SetLineColor(1)
    elif Sample=='QCD' or Sample=="QCDMuEnriched":
        hists[Sample].SetFillColor(ROOT.kCyan-1)
        hists[Sample].SetLineColor(1)
    else:
        hists[Sample].SetFillColor(ROOT.kGray+3)
        hists[Sample].SetLineColor(1)
    
    bg.Add(hists[Sample])
    l.AddEntry(hists[Sample], Sample.replace("NLO", ""), 'f')

for i, Sample in enumerate(Samples):
    print(i, Sample)
    if i == 0:
        h_bg = hists[Sample].Clone("Summed")
    else:
        h_bg.Add(hists[Sample])

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

xRangeLow = 0
xRangeHigh = 250
if whichVar == 'met':
    xRangeHigh = 350
elif whichVar == 'jetPt':
    xRangeLow = 20
    xRangeHigh = 600

padUp.Draw()
padUp.cd()
padUp.SetLogy()
data.Rebin(rebin)
#bg.Draw('hist')
bg.GetXaxis().SetRangeUser(xRangeLow,xRangeHigh)
bg.SetMinimum(100)
data.Draw()
data.SetMarkerSize(1)
data.SetMarkerStyle(20)
data.SetLineColor(1)
l.Draw('same')
#bg.SetMinimum(1000)
p1,p2,p3,p4,p5,p6,p7=RatioPlotTopXaxisSetup
SetupAxis(data, p1,p2,p3,p4,p5,p6,p7)
#h4.GetYaxis().SetMaxDigits(2)
p1,p2,p3,p4,p5,p6,p7=RatioPlotTopYaxisSetup
p5='N'
SetupAxis(data, p1,p2,p3,p4,p5,p6,p7)


c.cd()
padDown.Draw()
padDown.cd()
#padDown.cd()
ratio = data.Clone("Ratio")
#h_bg = bg.Sum()
ratio.Divide(h_bg)
ratio.Draw()
 
ratio.GetXaxis().SetRangeUser(xRangeLow,xRangeHigh)
#ratio.GetXaxis().SetTickSize(0.4)
ratio.SetMaximum(2)
ratio.SetMinimum(0)
p1,p2,p3,p4,p5,p6,p7=RatioPlotBottomXaxisSetup
if whichVar == 'mll':
    p5='m_{ll}'
elif whichVar == 'met':
    p5='MET'
elif whichVar == 'dRll':
    p5 = 'dR_{ll}'
elif whichVar == 'jetPt':
    p5 = 'p_{t,jet}'
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
c.SaveAs('data_mc_{}_{}_{}.pdf'.format(whichVar, whichTree, whichFolder))
        
            
