import os, sys
import ROOT

fil1 = ROOT.TFile("crabCfg/MuMuNtuple_SingleMuon_Run2016B.root")
fil2 = ROOT.TFile("crabCfg/MuMuNtuple_DoubleMuon_Run2016B.root")

events1 = fil1.Get("analyzeMuonsForEmbedding/analysisTree")
events2 = fil2.Get("analyzeMuonsForEmbedding/analysisTree")

out = ROOT.TFile("comparisonOut.root", "recreate")

h={}
h['dR1'] = ROOT.TH1F("dR1", "dR", 65, 0, 6.5)
h['dR1_jet500'] = ROOT.TH1F("dR1_jet500", "dR", 65, 0, 6.5)
h['dR2'] = ROOT.TH1F("dR2", "dR2", 65, 0, 6.5)
h['dR2_jet500'] = ROOT.TH1F("dR2_jet500", "dR2", 65, 0, 6.5)
h['m1'] = ROOT.TH1F("m1", "m1", 300, 0, 150)
h['m1_jet500'] = ROOT.TH1F("m1_jet500", "dR", 300, 0, 150)
h['m2'] = ROOT.TH1F("m2", "m2", 300, 0, 150)
h['m2_jet500'] = ROOT.TH1F("m2_jet500", "m2", 300, 0, 150)


N = events1.GetEntries()
for i in range(N):
    events1.GetEntry(i)
    dR = getattr(events1, 'mumuInfoBoosted/dR')
    mass = getattr(events1, 'mumuInfoBoosted/mass')
    h['dR1'].Fill(dR)
    h['m1'].Fill(mass)
    
    jetPt = getattr(events1, 'jetInfo/pt')
    if jetPt > 500:
        dR = getattr(events1, 'mumuInfoBoosted/dR')
        mass = getattr(events1, 'mumuInfoBoosted/mass')
        h['dR1_jet500'].Fill(dR)
        h['m1_jet500'].Fill(mass)

N = events2.GetEntries()
for i in range(N):
    events2.GetEntry(i)
    dR = getattr(events2, 'mumuInfoBoosted/dR')
    mass = getattr(events2, 'mumuInfoBoosted/mass')
    h['dR2'].Fill(dR)
    h['m2'].Fill(mass)
    jetPt = getattr(events2, 'jetInfo/pt')
    if jetPt > 500:
        dR = getattr(events2, 'mumuInfoBoosted/dR')
        mass = getattr(events2, 'mumuInfoBoosted/mass')
        h['dR2_jet500'].Fill(dR)
        h['m2_jet500'].Fill(mass)

#c = ROOT.TCanvas("Comparisons", "Comparisons", 0, 0, 650, 800)


out.cd()
for key in h.keys():
    h[key].Write()
out.Close()
