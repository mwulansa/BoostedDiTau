import ROOT, sys, os
import numpy as np

ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/JetInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/MuonInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/ElectronInfoDS.h"')

inputFileListName=sys.argv[1]
inputFileList=inputFileListName

if len(sys.argv)>2:
    outputFileDir=sys.argv[2]
else:
    outputFileDir = "./output/"

#outputFileDir = "./output/"
outputFileName = outputFileDir+"h_Baseline_EMu_"+inputFileListName.split("/")[-1].replace(".txt",".root")

#out=ROOT.TFile.Open("testEMu.root",'recreate')
out=ROOT.TFile.Open(outputFileName,'recreate')

fchain = ROOT.TChain('tcpNtuples/analysisTree')

#gchain = ROOT.TChain('tcpGenNtuples/genTree')

pi = np.pi

h = {}

h['hEvents'] = ROOT.TH1F ("NEvents", "Number of Events; ;N", 2, 0, 2)

h['hJetPt'] = ROOT.TH1F ("hJetPt", "Jet P_{T} ; P_{T} ; N", 2000, 0, 2000)
h['hBtag'] = ROOT.TH1F ("hBtag", "DeepJet Score; score; N", 50, 0, 1)
h['hBJetPt'] = ROOT.TH1F ("hBJetPt", "B-tagged Jet P_{T} ; P_{T} ; N", 2000, 0, 2000)
h['hMuonPt'] = ROOT.TH1F ("hMuPt", "Muon P_{T} ; P_{T} ; N", 500, 0, 500)
h['muonId'] = ROOT.TH1F ("muonId", ";;", 4, 0, 4)
h['hElectronPt'] = ROOT.TH1F ("hEPt", "Electron P_{T} ; P_{T} ; N", 500, 0, 500)
h['hMetPt'] = ROOT.TH1F ("hMetPt", "MET P_{T} ; P_{T} ; N", 500, 0, 500)

h['hEMuBaseline'] = ROOT.TH1F ("hEMu_Baseline", "Baseline Visible Mass;M_{e#mu};N", 500, 0, 100)
h['hEMudRcut'] = ROOT.TH1F ("hEMu_dRcut", "After dR cut Visible Mass;M_{e#mu};N", 500, 0, 100)
h['hEMuMetcut'] = ROOT.TH1F ("hEMu_Metcut", "After MetPt cut Visible Mass;M_{e#mu};N", 500, 0, 100)

h['hEMuBaseline_dR'] = ROOT.TH1F ("hEMu_Baseline_dR", "Baseline DeltaR(e,#mu);dR_{e,#mu};N", 100, 0, 5)
h['hEMuBaseline_dRj'] = ROOT.TH1F ("hEMu_Baseline_dRj", "Baseline DeltaR(l,jet);dR_{l,jet};N", 100, 0, 5)
h['hEMuBaseline_MetPt'] = ROOT.TH1F ("hEMu_Baseline_MetPt", "Baseline Met P_{T}; P_{T} ;N", 500, 0, 500)
h['hEMuBaseline_dPhi'] = ROOT.TH2F ("hEMu_Baseline_dPhi", "Baseline dPhi lepton-jet; dPhi(E_{T}, #mu); dPhi(E_{T}, jet)", 100, -pi, pi, 100, -pi, pi)
h['hEMuBaseline_dPhil'] = ROOT.TH2F ("hEMu_Baseline_dPhil", "Baseline dPhi leptons; dPhi(E_{T}, #mu); dPhi(E_{T}, e)", 100, -pi, pi, 100, -pi, pi)

h['hEMuSelection'] = ROOT.TH1F ("hEMu_AfterSelection", "Mass;dR_{e#mu};N Events", 500, 0, 100)

for key in h.keys():
    h[key].Sumw2()

#prefix = "root://cmseos.fnal.gov//store/user/mwulansa/TCPNtuple/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/Ntuple_DYJetsToLL_M-50_Summer20UL17_v5-2/211117_035810/0000/"

#prefix = "root://cmseos.fnal.gov//store/user/mwulansa/TCPNtuple/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/Ntuple_DYJetsToLL_M-50_Summer20UL17_v5-2/211117_035810/0000/"

inputFileNames=open(inputFileList, 'r')
for inputFileName in inputFileNames:
    inputFileName=inputFileName.replace("\n","")
    print(inputFileName.replace("\n",""))

    fchain.Add(inputFileName)
    fchain.AddFriend('tcpTrigNtuples/triggerTree', inputFileName)
    fchain.AddFriend('lumiSummary/lumiTree', inputFileName)

#    print(fchain.GetEntries())
#    gchain.Add(inputFileName)

    jets = ROOT.JetInfoDS()
    muons = ROOT.MuonInfoDS()
    electrons = ROOT.ElectronInfoDS()

    fchain.SetBranchAddress("Jets", ROOT.AddressOf(jets))
    fchain.SetBranchAddress("Muons", ROOT.AddressOf(muons))
    fchain.SetBranchAddress("Electrons", ROOT.AddressOf(electrons))

    for iev in range(fchain.GetEntries()): # Be careful!!!                                                                                                   
       fchain.GetEntry(iev)

       mets = fchain.GetBranch("Mets")
       met_pt = mets.GetLeaf('pt').GetValue()
       met_phi = mets.GetLeaf('phi').GetValue()

       weight = fchain.GetBranch("lumiInfo")
       genweight = weight.GetLeaf('weight').GetValue()

       h['hEvents'].Fill(0.5, 1)
       h['hEvents'].Fill(1.5, genweight)

       h['hMetPt'].Fill(met_pt)

       isSingleJet = fchain.GetLeaf('isSingleJet').GetValue()
       isHT = fchain.GetLeaf('isHT').GetValue()
       isHTMHT = fchain.GetLeaf('isHTMHT').GetValue()
       isMu = fchain.GetLeaf('isMu').GetValue()
       isIsoMu = fchain.GetLeaf('isIsoMu').GetValue()
       isIsoMuTau = fchain.GetLeaf('isIsoMuTau').GetValue()
       isIsoEle = fchain.GetLeaf('isIsoEle').GetValue()
       isEleTau = fchain.GetLeaf('isEleTau').GetValue()
       isMuonEG = fchain.GetLeaf('isMuonEG').GetValue()

       s_j = []
       s_b = []
       s_e = []
       s_mu = []

       if jets.size()>0:
          for i in range(jets.size()):
             jet = jets.at(i)
             if jet.id == 2:
                 h['hJetPt'].Fill(jet.pt)
                 h['hBtag'].Fill(jet.deepjet)
                 s_j+=[jet]
                 if jet.deepjet > 0.7476:
                     h['hBJetPt'].Fill(jet.pt)
                     s_b+=[jet]

       if muons.size()>0:
          for i in range(muons.size()):
             muon = muons.at(i)
             if muon.id >= 1:
                 h['hMuonPt'].Fill(muon.pt) 
                 s_mu+=[muon]

       if electrons.size()>0:
          for i in range(electrons.size()):
             electron = electrons.at(i)
             if electron.id >= 1 :
                 h['hElectronPt'].Fill(electron.pt)
                 s_e+=[electron]

       s_j.sort(key=lambda x: x.pt, reverse=True)
       s_e.sort(key=lambda x: x.pt, reverse=True)
       s_mu.sort(key=lambda x: x.pt, reverse=True)

#       if isHT == 1 or isSingleJet == 1 or isMu == 1 or \\
#          isIsoMu == 1 or isIsoMuTau == 1 or \\
#          isIsoEle == 1 or isEleTau == 1 or \\
#          isMuonEG == 1: 

       if len(s_e) > 0 and len(s_j) > 0 and len(s_mu) > 0 and s_e[0].charge*s_mu[0].charge < 0 and len(s_b) == 0:

           mu = ROOT.TLorentzVector()
           mu.SetPtEtaPhiM(s_mu[0].pt, s_mu[0].eta, s_mu[0].phi, s_mu[0].mass)

           e = ROOT.TLorentzVector()
           e.SetPtEtaPhiM(s_e[0].pt, s_e[0].eta, s_e[0].phi, s_e[0].mass)

           j = ROOT.TLorentzVector()
           j.SetPtEtaPhiM(s_j[0].pt, s_j[0].eta, s_j[0].phi, s_j[0].mass)

           m = ROOT.TLorentzVector()
           m.SetPtEtaPhiM(met_pt, 0, met_phi, 0)

           h['hEMuBaseline'].Fill((e+mu).M(), genweight)

           if mu.DeltaR(e) < 0.4 and mu.DeltaR(j) > 0.8 and e.DeltaR(j) > 0.8:
               h['hEMudRcut'].Fill((e+mu).M(), genweight)
               if met_pt > 100:
                   h['hEMuMetcut'].Fill((e+mu).M(), genweight)
                   if abs(m.DeltaPhi(mu))<1 and abs(m.DeltaPhi(j))>2:
                       h['hEMuSelection'].Fill((e+mu).M(), genweight)

           if met_pt > 100:
               if abs(m.DeltaPhi(mu))<1 and abs(m.DeltaPhi(j))>2:
                   h['hEMuBaseline_dR'].Fill(mu.DeltaR(e), genweight)
                   h['hEMuBaseline_dRj'].Fill((e+mu).DeltaR(j), genweight)

           if mu.DeltaR(e) < 0.4 and mu.DeltaR(j) > 0.8 and e.DeltaR(j) > 0.8:
               if met_pt > 100:
                   h['hEMuBaseline_dPhi'].Fill(m.DeltaPhi(mu), m.DeltaPhi(j), genweight)
                   h['hEMuBaseline_dPhil'].Fill(m.DeltaPhi(mu), m.DeltaPhi(e), genweight)

           if mu.DeltaR(e) < 0.4 and mu.DeltaR(j) > 0.8 and e.DeltaR(j) > 0.8:
               if abs(m.DeltaPhi(mu))<1 and abs(m.DeltaPhi(j))>2:
                   h['hEMuBaseline_MetPt'].Fill(met_pt, genweight)


out.cd()

for key in h.keys():
    h[key].Write()

out.Close()
