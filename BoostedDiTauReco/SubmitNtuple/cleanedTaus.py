import ROOT, sys, os
import numpy as np

ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/JetInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/MuonInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/ElectronInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/TauInfoDS.h"')

inputFileListName=sys.argv[1]
inputFileList=inputFileListName

#out=ROOT.TFile.Open("testEMu.root",'recreate')
out=ROOT.TFile.Open("studyTauClean.root",'recreate')

fchain = ROOT.TChain('tcpNtuples/analysisTree')

#gchain = ROOT.TChain('tcpGenNtuples/genTree')

pi = np.pi

h = {}

h['hEvents'] = ROOT.TH1F ("NEvents", "Number of Events; ;N", 2, 0, 2)

h['hJetPt'] = ROOT.TH1F ("hJetPt", "Jet P_{T} ; P_{T} ; N", 2000, 0, 2000)
h['hBtag'] = ROOT.TH1F ("hBtag", "DeepJet Score; score; N", 50, 0, 1)
h['hBJetPt'] = ROOT.TH1F ("hBJetPt", "B-tagged Jet P_{T} ; P_{T} ; N", 2000, 0, 2000)
h['hMuonPt'] = ROOT.TH1F ("hMuPt", "Muon P_{T} ; P_{T} ; N", 500, 0, 500)
h['hElectronPt'] = ROOT.TH1F ("hEPt", "Electron P_{T} ; P_{T} ; N", 500, 0, 500)
h['hMetPt'] = ROOT.TH1F ("hMetPt", "MET P_{T} ; P_{T} ; N", 500, 0, 500)

h['hTauUnCleaned'] = ROOT.TH1F ("hTauUnCleaned", "Taus_Uncleaned P_{T}; P_{T}; a.u.", 500, 0, 500)
h['hTauECleaned'] = ROOT.TH1F ("hTauECleaned", "Taus_ECleaned P_{T}; P_{T}; a.u.", 500, 0, 500)
h['hTauMuCleaned'] = ROOT.TH1F ("hTauMuCleaned", "Taus_MuCleaned P_{T}; P_{T}; a.u.", 500, 0, 500)
h['hTauLowCleaned'] = ROOT.TH1F ("hTauLowCleaned", "Taus_LowCleaned P_{T}; P_{T}; a.u.", 500, 0, 500)
h['hTauBoosted'] = ROOT.TH1F ("hTauBoosted", "Taus_Boosted P_{T}; P_{T}; a.u.", 500, 0, 500)

h['hEMuBaseline_munclean'] = ROOT.TH1F ("hEMu_Baseline_MUnCleaned", "Baseline Visible Mass;M_{#mu#tau};N", 100, 0, 100)
h['hEMuBaseline_eunclean'] = ROOT.TH1F ("hEMu_Baseline_EUnCleaned", "Baseline Visible Mass;M_{e#tau};N", 100, 0, 100)
h['hEMuBaseline_loweunclean'] = ROOT.TH1F ("hEMu_Baseline_LowPtEUnCleaned", "Baseline Visible Mass;M_{e#tau};N", 100, 0, 100)
h['hEMuBaseline_eclean'] = ROOT.TH1F ("hEMu_Baseline_ECleaned", "Baseline Visible Mass;M_{e#tau};N", 100, 0, 100)
h['hEMuBaseline_mclean'] = ROOT.TH1F ("hEMu_Baseline_MuCleaned", "Baseline Visible Mass;M_{#mu#tau};N", 100, 0, 100)
h['hEMuBaseline_lowclean'] = ROOT.TH1F ("hEMu_Baseline_LowPtECleaned", "Baseline Visible Mass;M_{e#tau};N", 100, 0, 100)

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
    tausUnCleaned = ROOT.TauInfoDS()
    tausECleaned = ROOT.TauInfoDS()
    tausMCleaned = ROOT.TauInfoDS()
    tausLowPtECleaned = ROOT.TauInfoDS()
    tausBoosted = ROOT.TauInfoDS()

    fchain.SetBranchAddress("Jets", ROOT.AddressOf(jets))
    fchain.SetBranchAddress("Muons", ROOT.AddressOf(muons))
    fchain.SetBranchAddress("Electrons", ROOT.AddressOf(electrons))
    fchain.SetBranchAddress("TausUnCleaned", ROOT.AddressOf(tausUnCleaned))
    fchain.SetBranchAddress("TausECleaned", ROOT.AddressOf(tausECleaned))
    fchain.SetBranchAddress("TausMCleaned", ROOT.AddressOf(tausMCleaned))
    fchain.SetBranchAddress("TausLowPtECleaned", ROOT.AddressOf(tausLowPtECleaned))
    fchain.SetBranchAddress("TausBoosted", ROOT.AddressOf(tausBoosted))

    print("jets", jets)
    print("taus", tausUnCleaned)

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

       unclean = []
       eclean = []
       mclean = []
       lowclean = []
       boosted = []

       if tausUnCleaned.size()>0:
           for i in range(tausUnCleaned.size()):
               tau = tausUnCleaned.at(i)
               if tau.mvaid >= 4:
                   unclean+=[tau]

       if tausECleaned.size()>0:
           for i in range(tausECleaned.size()):
               tau = tausECleaned.at(i)
               if tau.mvaid >= 4:
                   eclean+=[tau]

       if tausMCleaned.size()>0:
           for i in range(tausMCleaned.size()):
               tau = tausMCleaned.at(i)
               if tau.mvaid >= 4:
                   mclean+=[tau]

       if tausLowPtECleaned.size()>0:
           for i in range(tausLowPtECleaned.size()):
               tau = tausLowPtECleaned.at(i)
               if tau.mvaid >= 4:
                   lowclean+=[tau]

       if tausBoosted.size()>0:
           for i in range(tausBoosted.size()):
               tau = tausBoosted.at(i)
               if tau.mvaid >= 1:
                   boosted+=[tau]

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
       unclean.sort(key=lambda x: x.pt, reverse=True)
       eclean.sort(key=lambda x: x.pt, reverse=True)
       mclean.sort(key=lambda x: x.pt, reverse=True)
       lowclean.sort(key=lambda x: x.pt, reverse=True)
       boosted.sort(key=lambda x: x.pt, reverse=True)


       #Tau_mu #Tau_h

       if len(s_mu) > 0 and len(mclean) > 0 and len(s_j) > 0 and len(s_b) == 0 and mclean[0].charge*s_mu[0].charge < 0 :
           mu = ROOT.TLorentzVector()
           mu.SetPtEtaPhiM(s_mu[0].pt, s_mu[0].eta, s_mu[0].phi, s_mu[0].mass)

           tau = ROOT.TLorentzVector()
           tau.SetPtEtaPhiM(mclean[0].pt, mclean[0].eta, mclean[0].phi, mclean[0].mass)

           j = ROOT.TLorentzVector()
           j.SetPtEtaPhiM(s_j[0].pt, s_j[0].eta, s_j[0].phi, s_j[0].mass)



       if len(unclean) > 0 and len(s_mu) > 0 and len(s_j) and len(s_b) == 0 and unclean[0].charge*s_mu[0].charge < 0 :

           mu = ROOT.TLorentzVector()
           mu.SetPtEtaPhiM(s_mu[0].pt, s_mu[0].eta, s_mu[0].phi, s_mu[0].mass)

           tau = ROOT.TLorentzVector()
           tau.SetPtEtaPhiM(unclean[0].pt, unclean[0].eta, unclean[0].phi, unclean[0].mass)

           j = ROOT.TLorentzVector()
           j.SetPtEtaPhiM(s_j[0].pt, s_j[0].eta, s_j[0].phi, s_j[0].mass)

           if mu.DeltaR(tau) < 0.4 and mu.DeltaR(j) > 0.8 and tau.DeltaR(j) > 0.8:
               h['hEMuBaseline_munclean'].Fill((mu+tau).M(), genweight)

       if len(mclean) > 0 and len(s_mu) > 0 and len(s_j) and len(s_b) == 0 and unclean[0].charge*s_mu[0].charge < 0 :

           mu = ROOT.TLorentzVector()
           mu.SetPtEtaPhiM(s_mu[0].pt, s_mu[0].eta, s_mu[0].phi, s_mu[0].mass)

           tau = ROOT.TLorentzVector()
           tau.SetPtEtaPhiM(unclean[0].pt, unclean[0].eta, unclean[0].phi, unclean[0].mass)

           j = ROOT.TLorentzVector()
           j.SetPtEtaPhiM(s_j[0].pt, s_j[0].eta, s_j[0].phi, s_j[0].mass)

           if mu.DeltaR(tau) < 0.4 and mu.DeltaR(j) > 0.8 and tau.DeltaR(j) > 0.8:
               h['hEMuBaseline_mclean'].Fill((mu+tau).M(), genweight)


out.cd()

for key in h.keys():
    h[key].Write()

out.Close()
