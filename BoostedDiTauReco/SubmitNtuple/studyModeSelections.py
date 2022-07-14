import ROOT, sys, os
import numpy as np
import time

start_time = time.time()

ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/JetInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/MuonInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/ElectronInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/TauInfoDS.h"')

inputFileListName=sys.argv[1]
inputFileList=inputFileListName

if len(sys.argv)>2:
    outputFileDir=sys.argv[2]
else:
    outputFileDir = "./output/"

outputFileName = outputFileDir+"h_ModeSelectionStudy_"+inputFileListName.split("/")[-1].replace(".txt",".root")

out=ROOT.TFile.Open(outputFileName,'recreate')
print(outputFileName)

fchain = ROOT.TChain('tcpNtuples/analysisTree')

#gchain = ROOT.TChain('tcpGenNtuples/genTree')

pi = np.pi

h = {}

h['hEvents'] = ROOT.TH1F ("NEvents", "Number of Events; ;N", 2, 0, 2)

h['hTauEvents'] = ROOT.TH1F ("hTauEvents", "Number of each decay modes;;a.u.",7,0,7)

h['hTauTauTrigger_MetHT_nLepton'] = ROOT.TH1F ("hTauTauTrigger_MetHT_nLepton", "Number of leptons in TauTau events;;a.u.", 2, 0, 2)
h['hTauTauTrigger_MetHT_ElectrondR'] = ROOT.TH2F ("hTauTauTrigger_MetHT_ElectrondR", "Delta R between electrons and taus;#Delta R(e,tau1);#Delta R(e,tau2)", 100, 0, 5, 100, 0, 5)
h['hTauTauTrigger_MetHT_MuondR'] = ROOT.TH2F ("hTauTauTrigger_MetHT_MuondR", "Delta R between muonons and taus;#Delta R(mu,tau1);#Delta R(mu,tau2)", 100, 0, 5, 100, 0, 5)

h['hMuMu_DeltaRl'] = ROOT.TH1F ("hMuMu_DeltaRl", "Delta R between muons ;#Delta R(mu1,mu2);a.u.", 100, 0, 5)
h['hMuMu_DeltaRj'] = ROOT.TH2F ("hMuMu_DeltaRj", "Delta R between muons and jet ;#Delta R(mu1,j);#Delta R(mu2,j)", 100, 0, 5, 100, 0, 5)

h['hJetPt'] = ROOT.TH1F ("hJetPt", "Jet P_{T} ; P_{T} ; N", 1000, 0, 1000)
h['hBtag'] = ROOT.TH1F ("hBtag", "DeepJet Score; score; N", 50, 0, 1)
h['hBJetPt'] = ROOT.TH1F ("hBJetPt", "B-tagged Jet P_{T} ; P_{T} ; N", 1000, 0, 1000)
h['hMuonPt'] = ROOT.TH1F ("hMuPt", "Muon P_{T} ; P_{T} ; N", 500, 0, 500)
h['hIsoMuonPt'] = ROOT.TH1F ("hIsoMuPt", "Isolated Muon P_{T} ; P_{T} ; N", 500, 0, 500)
h['hElectronPt'] = ROOT.TH1F ("hEPt", "Electron P_{T} ; P_{T} ; N", 500, 0, 500)
h['hIsoElectronPt'] = ROOT.TH1F ("hIsoEPt", "Isolated Electron P_{T} ; P_{T} ; N", 500, 0, 500)
h['hMetPt'] = ROOT.TH1F ("hMetPt", "MET P_{T} ; P_{T} ; N", 500, 0, 500)

h['hTauUnCleanedPt'] = ROOT.TH1F ("hTauUnCleanedPt", "Uncleaned Tau P_{T}; P_{T}; a.u.", 500, 0, 500)
h['hMuCleanedPt'] = ROOT.TH1F ("hMuCleanedPt", "Mu-cleaned Tau P_{T}; P_{T}; a.u.", 500, 0, 500)
h['hECleanedPt'] = ROOT.TH1F ("hECleanedPt", "E-cleaned Tau P_{T}; P_{T}; a.u.", 500, 0, 500)
h['hTauBoostedPt'] = ROOT.TH1F ("hTauBoostedPt", "Boosted Tau P_{T}; P_{T}; a.u.", 500, 0, 500)

h['hTauUnCleaned'] = ROOT.TH1F ("hTauUnCleaned", "Taus_Uncleaned P_{T}; P_{T}; a.u.", 500, 0, 500)
h['hTauECleaned'] = ROOT.TH1F ("hTauECleaned", "Taus_ECleaned P_{T}; P_{T}; a.u.", 500, 0, 500)
h['hTauMuCleaned'] = ROOT.TH1F ("hTauMuCleaned", "Taus_MuCleaned P_{T}; P_{T}; a.u.", 500, 0, 500)
#h['hTauLowCleaned'] = ROOT.TH1F ("hTauLowCleaned", "Taus_LowCleaned P_{T}; P_{T}; a.u.", 500, 0, 500)
h['hTauBoosted'] = ROOT.TH1F ("hTauBoosted", "Taus_Boosted P_{T}; P_{T}; a.u.", 500, 0, 500)

h['hEMuBaseline'] = ROOT.TH1F ("hEMu_Baseline", "EMu Baseline ;M_{#tau#tau};a.u.", 100, 0, 100)
h['hMuMuBaseline'] = ROOT.TH1F ("hMuMu_Baseline", "MuMu Baseline ;M_{#tau#tau};a.u.", 100, 0, 100)
h['hMuTauBaseline'] = ROOT.TH1F ("hMuTau_Baseline", "MuTau Baseline ;M_{#tau#tau};a.u.", 100, 0, 100)
h['hETauBaseline'] = ROOT.TH1F ("hETau_Baseline", "ETau Baseline ;M_{#tau#tau};a.u.", 100, 0, 100)
h['hTauTauBaseline'] = ROOT.TH1F ("hTauTau_Baseline", "TauTau Baseline ;M_{#tau#tau};a.u.", 60, 0, 100)

h['hEESelection'] = ROOT.TH1F ("hEE_Selection", "EE After Selection ;M_{#tau#tau};a.u.", 100, 0, 100)
h['hEMuSelection'] = ROOT.TH1F ("hEMu_Selection", "EMu After Selection ;M_{#tau#tau};a.u.", 100, 0, 100)
h['hMuMuSelection'] = ROOT.TH1F ("hMuMu_Selection", "MuMu After Selection ;M_{#tau#tau};a.u.", 100, 0, 100)
h['hMuTauSelection'] = ROOT.TH1F ("hMuTau_Selection", "MuTau After Selection ;M_{#tau#tau};a.u.", 100, 0, 100)
h['hETauSelection'] = ROOT.TH1F ("hETau_Selection", "ETau After Selection ;M_{#tau#tau};a.u.", 100, 0, 100)
h['hTauTauSelection'] = ROOT.TH1F ("hTauTau_Selection", "TauTau After Selection ;M_{#tau#tau};a.u.", 100, 0, 100)

#EE-------------------------------

h['hEE_JetPt'] = ROOT.TH1F ("hEE_JetPt", "Jet P_{T} ; P_{T} ; a.u.", 1000, 0, 1000)
h['hEE_MetPt'] = ROOT.TH1F ("hEE_MetPt", "Met P_{T} ; P_{T} ; a.u.", 500, 0, 500)
h['hEE_ElectronPt'] = ROOT.TH1F ("hEE_ElectronPt", "Muon P_{T} ; P_{T} ; a.u.", 500, 0, 500)

h['hEESelection_JetPt'] = ROOT.TH1F ("hEE_Selection_JetPt", "Jet P_{T} ; P_{T} ; a.u.", 1000, 0, 1000)
h['hEESelection_MetPt'] = ROOT.TH1F ("hEE_Selection_MetPt", "Met P_{T} ; P_{T} ; a.u.", 500, 0, 500)
h['hEESelection_ElectronPt'] = ROOT.TH1F ("hEE_Selection_ElectronPt", "Electrion P_{T} ; P_{T} ; a.u.", 500, 0, 500)

h['hEEBaseline_dRlj'] = ROOT.TH2F ("hEE_Baseline_dRlj", "Baseline DeltaR lepton/jet ; dR_{e1,jet} ; dR_{e2,jet}", 100, 0, 5, 100, 0, 5)
h['hEEBaseline_dRl'] = ROOT.TH2F ("hEE_Baseline_dRl", "Baseline DeltaR(e1,e2) ; M_{#tau#tau} ; dR_{e1,e2}", 100, 0, 100, 100, 0, 5)
h['hEEBaseline_dRj'] = ROOT.TH2F ("hEE_Baseline_dRj", "Baseline DeltaR(lep,jet) ; M_{#tau#tau} ; dR_{e1,jet}", 100, 0, 100, 100, 0, 5)
h['hEEBaseline_dR'] = ROOT.TH2F ("hEE_Baseline_dR", "Baseline DeltaR ; dR_{e1,e2} ; dR_{e1,jet}", 100, 0, 5, 100, 0, 5)

h['hEEdRcut_dRl'] = ROOT.TH2F ("hEE_dRcut_dRl", "After dR DeltaR(e1,e2) ; M_{#tau#tau} ; dR_{e1,e2}", 100, 0, 100, 100, 0, 5)
h['hEEdRcut_dRj'] = ROOT.TH2F ("hEE_dRcut_dRj", "After dR DeltaR(e1,e2) ; M_{#tau#tau} ; dR_{e1,jet}", 100, 0, 100, 100, 0, 5)

h['hEEBaseline_dPhij'] = ROOT.TH2F ("hEE_Baseline_dPhi", "Baseline dPhi lepton/jet; dPhi(E_{T}, #mu); dPhi(E_{T}, jet)", 100, -pi, pi, 100, -pi, pi)
h['hEEBaseline_dPhil'] = ROOT.TH2F ("hEE_Baseline_dPhil", "Baseline dPhi leptons; dPhi(E_{T}, #mu); dPhi(E_{T}, e)", 100, -pi, pi, 100, -pi, pi)

h['hEETrigger_dRl'] = ROOT.TH2F ("hEE_Trigger_dRl", "After Trigger DeltaR leptons ; M_{#tau#tau} ; dR_{e1,e2}", 100, 0, 100, 100, 0, 5)
h['hEETrigger_dRj'] = ROOT.TH2F ("hEE_Trigger_dRj", "After Trigger DeltaR leptons/jet ; M_{#tau#tau} ; dR_{l,j}", 100, 0, 100, 100, 0, 5)
h['hEETrigger_dPhil'] = ROOT.TH2F ("hEE_Trigger_dPhil", "After Trigger dPhi lepton; M_{#tau#tau}; dPhi(E_{T}, l)", 100, 0, 100, 100, -pi, pi)
h['hEETrigger_dPhij'] = ROOT.TH2F ("hEE_Trigger_dPhij", "After Trigger dPhi jet; M_{#tau#tau}; dPhi(E_{T}, jet)", 100, 0, 100, 100, -pi, pi)
h['hEETrigger_MetPt'] = ROOT.TH2F ("hEE_Trigger_MetPt","After Trigger Met Pt;M_{#tau#tau};Met P_{T}", 100, 0, 100, 500, 0, 500)

#MuMu-----------------------------

h['hMuMu_JetPt'] = ROOT.TH1F ("hMuMu_JetPt", "Jet P_{T} ; P_{T} ; a.u.", 1000, 0, 1000)
h['hMuMu_MetPt'] = ROOT.TH1F ("hMuMu_MetPt", "Met P_{T} ; P_{T} ; a.u.", 500, 0, 500)
h['hMuMu_MuonPt'] = ROOT.TH1F ("hMuMu_MuonPt", "Muon P_{T} ; P_{T} ; a.u.", 500, 0, 500)

h['hMuMuSelection_JetPt'] = ROOT.TH1F ("hMuMu_Selection_JetPt", "Jet P_{T} ; P_{T} ; a.u.", 1000, 0, 1000)
h['hMuMuSelection_MetPt'] = ROOT.TH1F ("hMuMu_Selection_MetPt", "Met P_{T} ; P_{T} ; a.u.", 500, 0, 500)
h['hMuMuSelection_MuonPt'] = ROOT.TH1F ("hMuMu_Selection_MuonPt", "Muon P_{T} ; P_{T} ; a.u.", 500, 0, 500)

#EMu------------------------------

h['hEMu_JetPt'] = ROOT.TH1F ("hEMu_JetPt", "Jet P_{T} ; P_{T} ; a.u.", 1000, 0, 1000)
h['hEMu_MetPt'] = ROOT.TH1F ("hEMu_MetPt", "Met P_{T} ; P_{T} ; a.u.", 500, 0, 500)
h['hEMu_MuonPt'] = ROOT.TH1F ("hEMu_MuonPt", "Muon P_{T} ; P_{T} ; a.u.", 500, 0, 500)
h['hEMu_ElectronPt'] = ROOT.TH1F ("hEMu_ElectronPt", "Electron P_{T} ; P_{T} ; a.u.", 500, 0, 500)

h['hEMuSelection_JetPt'] = ROOT.TH1F ("hEMu_Selection_JetPt", "Jet P_{T} ; P_{T} ; a.u.", 1000, 0, 1000)
h['hEMuSelection_MetPt'] = ROOT.TH1F ("hEMu_Selection_MetPt", "Met P_{T} ; P_{T} ; a.u.", 500, 0, 500)
h['hEMuSelection_MuonPt'] = ROOT.TH1F ("hEMu_Selection_MuonPt", "Muon P_{T} ; P_{T} ; a.u.", 500, 0, 500)
h['hEMuSelection_ElectronPt'] = ROOT.TH1F ("hEMu_Selection_ElectronPt", "Electron P_{T} ; P_{T} ; a.u.", 500, 0, 500)

#MuTau--------------------------------

h['hMuTau_JetPt'] = ROOT.TH1F ("hMuTau_JetPt", "Jet P_{T} ; P_{T} ; a.u.", 1000, 0, 1000)
h['hMuTau_MetPt'] = ROOT.TH1F ("hMuTau_MetPt", "Met P_{T} ; P_{T} ; a.u.", 500, 0, 500)
h['hMuTau_MuonPt'] = ROOT.TH1F ("hMuTau_MuonPt", "Muon P_{T} ; P_{T} ; a.u.", 500, 0, 500)

h['hMuTauSelection_JetPt'] = ROOT.TH1F ("hMuTau_Selection_JetPt", "Jet P_{T} ; P_{T} ; a.u.", 1000, 0, 1000)
h['hMuTauSelection_MetPt'] = ROOT.TH1F ("hMuTau_Selection_MetPt", "Met P_{T} ; P_{T} ; a.u.", 500, 0, 500)
h['hMuTauSelection_MuonPt'] = ROOT.TH1F ("hMuTau_Selection_MuonPt", "Muon P_{T} ; P_{T} ; a.u.", 500, 0, 500)

h['hMuTauTrigger_dRl'] = ROOT.TH2F ("hMuTau_Trigger_dRl", "Trigger DeltaR(mu,tau) ; M_{#mu#tau} ; dR(mu,tau)", 100,0,100,100, 0, 0.5)

#ETau--------------------------------

h['hETau_JetPt'] = ROOT.TH1F ("hETau_JetPt", "Jet P_{T} ; P_{T} ; a.u.", 1000, 0, 1000)
h['hETau_MetPt'] = ROOT.TH1F ("hETau_MetPt", "Met P_{T} ; P_{T} ; a.u.", 500, 0, 500)
h['hETau_ElectronPt'] = ROOT.TH1F ("hETau_ElectronPt", "Electron P_{T} ; P_{T} ; a.u.", 500, 0, 500)

h['hETauSelection_JetPt'] = ROOT.TH1F ("hETau_Selection_JetPt", "Jet P_{T} ; P_{T} ; a.u.", 1000, 0, 1000)
h['hETauSelection_MetPt'] = ROOT.TH1F ("hETau_Selection_MetPt", "Met P_{T} ; P_{T} ; a.u.", 500, 0, 500)
h['hETauSelection_ElectronPt'] = ROOT.TH1F ("hETau_Selection_ElectronPt", "Electron P_{T} ; P_{T} ; a.u.", 500, 0, 500)

h['hETauTrigger_dRl'] = ROOT.TH2F ("hETau_Trigger_dRl", "Trigger DeltaR(e,tau) ; M_{e#tau} ; dR(e,tau)", 100,0,100,100, 0, 0.5)

h['hETauTrigger_JetPt'] = ROOT.TH1F ("hETau_Trigger_JetPt", "Jet P_{T} ; P_{T} ; a.u.", 1000, 0, 1000)
h['hETauTrigger_MetPt'] = ROOT.TH1F ("hETau_Trigger_MetPt", "Met P_{T} ; P_{T} ; a.u.", 500, 0, 500)

h['hETauTrigger_JetPt_MetPt'] = ROOT.TH2F ("hETau_Trigger_JetPt_MetPt", "JetPt_MetPt ; Jet P_{T} ; Met P_{T}", 1000, 0, 1000, 500, 0, 500)

#TauTau-------------------------

h['hTauTau_JetPt'] = ROOT.TH1F ("hTauTau_JetPt", "Jet P_{T} ; P_{T} ; a.u.", 1000, 0, 1000)
h['hTauTau_MetPt'] = ROOT.TH1F ("hTauTau_MetPt", "Met P_{T} ; P_{T} ; a.u.", 500, 0, 500)

h['hTauTauSelection_JetPt'] = ROOT.TH1F ("hTauTau_Selection_JetPt", "Jet P_{T} ; P_{T} ; a.u.", 1000, 0, 1000)
h['hTauTauSelection_MetPt'] = ROOT.TH1F ("hTauTau_Selection_MetPt", "Met P_{T} ; P_{T} ; a.u.", 500, 0, 500)

h['hTauTauTrigger_JetPt'] = ROOT.TH1F ("hTauTau_Trigger_JetPt", "Jet P_{T} ; P_{T} ; a.u.", 1000, 0, 1000)
h['hTauTauTrigger_MetPt'] = ROOT.TH1F ("hTauTau_Trigger_MetPt", "Met P_{T} ; P_{T} ; a.u.", 500, 0, 500)

h['hTauTauTrigger_JetPt_MetPt'] =ROOT.TH2F ("hTauTau_Trigger_JetPt_MetPt", "JetPt_MetPt ; Jet P_{T} ; Met P_{T}", 1000, 0, 1000, 500, 0, 500)

#Masses-----------------------------------

h['hEEBaseline'] = ROOT.TH1F ("hEE_Baseline", "EE Baseline ;M_{#tau#tau};a.u.", 100, 0, 100)
h['hEEdR'] = ROOT.TH1F ("hEE_dR", "EE After dR ;M_{#tau#tau};a.u.", 100, 0, 100)
h['hEEMetcut'] = ROOT.TH1F ("hEE_Metcut", "EE After Met cut ;M_{#tau#tau};a.u.", 100, 0, 100)

h['hEETrigger_SingleE'] = ROOT.TH1F ("hEE_Trigger_SingleE", "EE After Trigger (+SingleE) ;M_{#tau#tau};a.u.", 100, 0, 100)

h['hMuMuTrigger_SingleMu'] = ROOT.TH1F ("hMuMu_Trigger_SingleMu", "MuMu After Trigger (+SingleMu) ;M_{#tau#tau};a.u.", 100, 0, 100)

h['hEMuTrigger_SingleL'] = ROOT.TH1F ("hEMu_Trigger_SingleL", "EMu After Trigger (+SingleL) ;M_{#tau#tau};a.u.", 100, 0, 100)

h['hMuTauTrigger_SingleMu'] = ROOT.TH1F ("hMuTau_Trigger_SingleMu", "MuTau After Trigger (+SingleMu) ;M_{#tau#tau};a.u.", 100, 0, 100)

h['hETauTrigger_MetHT_Tight'] = ROOT.TH1F ("hETau_Trigger_MetHT_Tight", "ETau After Trigger (+MetHT) (Tight) ;M_{#tau#tau};a.u.", 100, 0, 100)
h['hETauTrigger_MetHT_Loose'] = ROOT.TH1F ("hETau_Trigger_MetHT_Loose", "ETau After Trigger (+MetHT) (Loose) ;M_{#tau#tau};a.u.", 100, 0, 100)
h['hETauTrigger_MetHT_Tight_wdR'] = ROOT.TH1F ("hETau_Trigger_MetHT_Tight_wdR", "ETau After Trigger (+MetHT) (Tight wdR) ;M_{#tau#tau};a.u.", 100, 0, 100)
h['hETauTrigger_MetHT_Loose_wdR'] = ROOT.TH1F ("hETau_Trigger_MetHT_Loose_wdR", "ETau After Trigger (+MetHT) (Loose wdR) ;M_{#tau#tau};a.u.", 100, 0, 100)

h['hTauTauTrigger_MetHT_Tight'] = ROOT.TH1F ("hTauTau_Trigger_MetHT_Tight", "TauTau After MetHT (Tight) ;M_{#tau#tau};a.u.", 100, 0, 100)
h['hTauTauTrigger_MetHT_Loose'] = ROOT.TH1F ("hTauTau_Trigger_MetHT_Loose", "TauTau After MetHT (loose) ;M_{#tau#tau};a.u.", 100, 0, 100)

h['hMuTauSelection_veto'] = ROOT.TH1F ("hMuTau_Selection_Veto", "MuTau After Selection - veto;M_{#tau#tau};a.u.", 100, 0, 100)
h['hETauSelection_veto'] = ROOT.TH1F ("hETau_Selection_Veto", "ETau After Selection - veto;M_{#tau#tau};a.u.", 100, 0, 100)
h['hTauTauSelection_veto'] = ROOT.TH1F ("hTauTau_Selection_Veto", "TauTau After Selection - veto;M_{#tau#tau};a.u.", 100, 0, 100)

h['hETauTrigger_EE'] = ROOT.TH1F ("hETau_Trigger_EE", "Boosted EE Events in ETau ;M_{#tau#tau};a.u.", 100, 0, 100)
h['hETauTrigger_EE_ETau'] = ROOT.TH1F ("hETau_Trigger_EE_ETau", "Compare EE and ETau ;M_{#tau#tau};a.u.", 100, 0, 100)

h['hTauTauTrigger_EE'] = ROOT.TH1F ("hTauTau_Trigger_EE", "Boosted EE Events in TauTau;M_{#tau#tau};a.u.", 100, 0, 100)
h['hTauTauTrigger_EE_TauTau'] = ROOT.TH1F ("hTauTau_Trigger_EE_TauTau", "Compare EE and TauTau ;M_{#tau#tau};a.u.", 100, 0, 100)
h['hTauTauTrigger_MuMu'] = ROOT.TH1F ("hTauTau_Trigger_MuMu", "Boosted MuMu Events in TauTau;M_{#tau#tau};a.u.", 100, 0, 100)
h['hTauTauTrigger_MuMu_TauTau'] = ROOT.TH1F ("hTauTau_Trigger_MuMu_TauTau", "Compare MuMu and TauTau ;M_{#tau#tau};a.u.", 100, 0, 100)


for key in h.keys():
    h[key].Sumw2()

#prefix = "root://cmseos.fnal.gov//store/user/mwulansa/TCPNtuple/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/Ntuple_DYJetsToLL_M-50_Summer20UL17_v5-2/211117_035810/0000/"

#prefix = "root://cmseos.fnal.gov//store/user/mwulansa/TCPNtuple/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/Ntuple_DYJetsToLL_M-50_Summer20UL17_v5-2/211117_035810/0000/"

def EE_Channel():

    isEE = 0 

    e1 = ROOT.TLorentzVector()
    e1.SetPtEtaPhiM(s_e[0].pt, s_e[0].eta, s_e[0].phi, s_e[0].mass)
    e2 = ROOT.TLorentzVector()
    e2.SetPtEtaPhiM(s_e[1].pt, s_e[1].eta, s_e[1].phi, s_e[1].mass)

    j = ROOT.TLorentzVector()
    j.SetPtEtaPhiM(s_j[0].pt, s_j[0].eta, s_j[0].phi, s_j[0].mass)

    m = ROOT.TLorentzVector()
    m.SetPtEtaPhiM(met_pt, 0, met_phi, 0)

    if (e1+e2).M() > 1:

        h['hEEBaseline'].Fill((e1+e2).M(), genweight)

        h['hEE_JetPt'].Fill(j.Pt(), genweight)
        h['hEE_MetPt'].Fill(met_pt, genweight)
        h['hEE_ElectronPt'].Fill(e1.Pt(), genweight)

        # h['hEEBaseline_dRl'].Fill((e1+e2).M(), e1.DeltaR(e2), genweight)
        # h['hEEBaseline_dRj'].Fill((e1+e2).M(), j.DeltaR(e1+e2), genweight)
        # h['hEEBaseline_dRlj'].Fill(e1.DeltaR(j), e2.DeltaR(j), genweight)

        # h['hEEBaseline_dR'].Fill(e1.DeltaR(e2), e1.DeltaR(j), genweight)

        # h['hEEBaseline_dPhil'].Fill(m.DeltaPhi(e1), m.DeltaPhi(e2), genweight)
        # h['hEEBaseline_dPhij'].Fill(m.DeltaPhi(e1), m.DeltaPhi(j), genweight)

        if e1.DeltaR(e2) < 0.4 and e1.DeltaR(j)> 0.8  and e2.DeltaR(j) > 0.8:
            h['hEEdR'].Fill((e1+e2).M(), genweight)

            h['hEEdRcut_dRl'].Fill((e1+e2).M(), e1.DeltaR(e2), genweight)
            h['hEEdRcut_dRj'].Fill((e1+e2).M(), j.DeltaR(e1+e2), genweight)

            if met_pt > 100:
                h['hEEMetcut'].Fill((e1+e2).M(), genweight)
                if abs(m.DeltaPhi(e1)) < 1 and abs(m.DeltaPhi(j)) > 2:
                    h['hEESelection'].Fill((e1+e2).M(), genweight)
                    h['hEESelection_JetPt'].Fill(j.Pt(), genweight)
                    h['hEESelection_MetPt'].Fill(met_pt, genweight)
                    h['hEESelection_ElectronPt'].Fill(e1.Pt(), genweight)

                    if ( j.Pt() > 500 and ( isSingleJet == 1 or isHT == 1 ) ) \
                       or ( j.Pt() > 500 and met_pt > 200 and isHTMHT == 1 ) \
                       or ( e1.Pt() > 35 and isIsoEle == 1 ):
                        h['hEETrigger_SingleE'].Fill((e1+e2).M(),  genweight)
                        # h['hEETrigger_dRl'].Fill((e1+e2).M(), e1.DeltaR(e2), genweight)
                        # h['hEETrigger_dRj'].Fill((e1+e2).M(), j.DeltaR(e1+e2), genweight)
                        # h['hEETrigger_MetPt'].Fill((e1+e2).M(), met_pt, genweight)
                        # h['hEETrigger_dPhil'].Fill((e1+e2).M(), m.DeltaPhi(e1), genweight)
                        # h['hEETrigger_dPhij'].Fill((e1+e2).M(), m.DeltaPhi(j), genweight)
                        isEE = 1
                        h['hTauEvents'].Fill(1,1)

    return isEE

def MuMu_Channel():

   isMuMu = 0

   mu1 = ROOT.TLorentzVector()
   mu1.SetPtEtaPhiM(s_mu[0].pt, s_mu[0].eta, s_mu[0].phi, s_mu[0].mass)

   mu2 = ROOT.TLorentzVector()
   mu2.SetPtEtaPhiM(s_mu[1].pt, s_mu[1].eta, s_mu[1].phi, s_mu[1].mass)

   j = ROOT.TLorentzVector()
   j.SetPtEtaPhiM(s_j[0].pt, s_j[0].eta, s_j[0].phi, s_j[0].mass)

   m = ROOT.TLorentzVector()
   m.SetPtEtaPhiM(met_pt, 0, met_phi, 0)

   if (mu1+mu2).M() > 1:

       h['hMuMuBaseline'].Fill((mu1+mu2).M(), genweight)

       h['hMuMu_JetPt'].Fill(j.Pt(), genweight)
       h['hMuMu_MetPt'].Fill(met_pt, genweight)
       h['hMuMu_MuonPt'].Fill(mu1.Pt(), genweight)

       h['hMuMu_DeltaRl'].Fill(mu1.DeltaR(mu2), genweight)
       h['hMuMu_DeltaRj'].Fill(mu1.DeltaR(j), mu2.DeltaR(j), genweight)

       if mu1.DeltaR(mu2) < 0.4 and mu1.DeltaR(j)> 0.8  and mu2.DeltaR(j) > 0.8:
           if met_pt > 100:
               if abs(m.DeltaPhi(mu1)) < 1 and abs(m.DeltaPhi(j)) > 2:

                   h['hMuMuSelection'].Fill((mu1+mu2).M(), genweight)
                   h['hMuMuSelection_JetPt'].Fill(j.Pt(), genweight)
                   h['hMuMuSelection_MetPt'].Fill(met_pt, genweight)
                   h['hMuMuSelection_MuonPt'].Fill(mu1.Pt(), genweight)

                   if ( j.Pt() > 500 and ( isSingleJet == 1 or isHT == 1 ) ) \
                      or ( mu1.Pt() > 50 and isMu == 1 ) or ( mu1.Pt() > 27 and isIsoMu == 1 ) :
                       h['hMuMuTrigger_SingleMu'].Fill((mu1+mu2).M(),  genweight)
                       isMuMu = 1
                       h['hTauEvents'].Fill(2,1)

   return isMuMu

def EMu_Channel():

   isEMu = 0

   mu = ROOT.TLorentzVector()
   mu.SetPtEtaPhiM(s_mu[0].pt, s_mu[0].eta, s_mu[0].phi, s_mu[0].mass)

   e = ROOT.TLorentzVector()
   e.SetPtEtaPhiM(s_e[0].pt, s_e[0].eta, s_e[0].phi, s_e[0].mass)

   j = ROOT.TLorentzVector()
   j.SetPtEtaPhiM(s_j[0].pt, s_j[0].eta, s_j[0].phi, s_j[0].mass)

   m = ROOT.TLorentzVector()
   m.SetPtEtaPhiM(met_pt, 0, met_phi, 0)

   if (mu+e).M() > 1:

       h['hEMuBaseline'].Fill((e+mu).M(), genweight)

       h['hEMu_JetPt'].Fill(j.Pt(), genweight)
       h['hEMu_MetPt'].Fill(met_pt, genweight)
       h['hEMu_MuonPt'].Fill(mu.Pt(), genweight)
       h['hEMu_ElectronPt'].Fill(e.Pt(), genweight)

       if mu.DeltaR(e) < 0.4 and mu.DeltaR(j)> 0.8  and e.DeltaR(j) > 0.8:
           if met_pt > 100:
               if abs(m.DeltaPhi(mu)) < 1 and abs(m.DeltaPhi(j)) > 2:

                   h['hEMuSelection'].Fill((mu+e).M(), genweight)
                   h['hEMuSelection_JetPt'].Fill(j.Pt(), genweight)
                   h['hEMuSelection_MetPt'].Fill(met_pt, genweight)
                   h['hEMuSelection_MuonPt'].Fill(mu.Pt(), genweight)
                   h['hEMuSelection_ElectronPt'].Fill(e.Pt(), genweight)

                   if ( j.Pt() > 500 and ( isSingleJet == 1 or isHT == 1 ) ) \
                      or ( isMuonEG == 1 and ( ( mu.Pt() > 8 and e.Pt() > 23 ) or ( mu.Pt() > 23 and e.Pt() > 12 ) ) ) \
                      or ( mu.Pt() > 50 and isMu == 1 ) or ( mu.Pt() > 27 and isIsoMu == 1 ) \
                      or ( e.Pt() > 35 and (isIsoEle == 1) ) :
                       h['hEMuTrigger_SingleL'].Fill((e+mu).M(),  genweight)
                       isEMu = 1
                       h['hTauEvents'].Fill(3,1)

   return isEMu

def MuTau_Channel():

       isMuTau = 0

       mu = ROOT.TLorentzVector()
       mu.SetPtEtaPhiM(s_mu[0].pt, s_mu[0].eta, s_mu[0].phi, s_mu[0].mass)

       tau = ROOT.TLorentzVector()
       tau.SetPtEtaPhiM(mclean[0].pt, mclean[0].eta, mclean[0].phi, mclean[0].mass)

       if len(s_j) > 1 :
           j0 = ROOT.TLorentzVector()
           j0.SetPtEtaPhiM(s_j[0].pt, s_j[0].eta, s_j[0].phi, s_j[0].mass)
           if tau.DeltaR(j0) < 0.3 :
               j = ROOT.TLorentzVector()
               j.SetPtEtaPhiM(s_j[1].pt, s_j[1].eta, s_j[1].phi, s_j[1].mass)
           else:
               j = j0
       if len(s_j) == 1:
           j = ROOT.TLorentzVector()
           j.SetPtEtaPhiM(s_j[0].pt, s_j[0].eta, s_j[0].phi, s_j[0].mass)

       m = ROOT.TLorentzVector()
       m.SetPtEtaPhiM(met_pt, 0, met_phi, 0)

       if (mu+tau).M() > 1:

           h['hMuTauBaseline'].Fill((mu+tau).M(), genweight)

           h['hMuTau_JetPt'].Fill(j.Pt(), genweight)
           h['hMuTau_MetPt'].Fill(met_pt, genweight)
           h['hMuTau_MuonPt'].Fill(mu.Pt(), genweight)

           if mu.DeltaR(tau) < 0.4 and mu.DeltaR(j)> 0.8  and tau.DeltaR(j) > 0.8:
               if met_pt > 100:
                   if abs(m.DeltaPhi(mu)) < 1 and abs(m.DeltaPhi(j)) > 2:

                       if len(s_isoe) == 0:
                           h['hMuTauSelection_veto'].Fill((mu+tau).M(), genweight)

                       h['hMuTauSelection'].Fill((mu+tau).M(), genweight)
                       h['hMuTauSelection_JetPt'].Fill(j.Pt(), genweight)
                       h['hMuTauSelection_MetPt'].Fill(met_pt, genweight)
                       h['hMuTauSelection_MuonPt'].Fill(mu.Pt(), genweight)

                       if ( j.Pt() > 500 and ( isSingleJet == 1 or isHT == 1) ) \
                          or ( ( mu.Pt() > 50 and isMu == 1 ) or ( mu.Pt() > 27 and isIsoMu == 1 ) ) :
                           h['hMuTauTrigger_SingleMu'].Fill((mu+tau).M(), genweight)
                           isMuTau = 1
                           h['hTauEvents'].Fill(4,1)

       return isMuTau

def ETau_Channel():

   isETauTight = 0

   e = ROOT.TLorentzVector()
   e.SetPtEtaPhiM(s_e[0].pt, s_e[0].eta, s_e[0].phi, s_e[0].mass)

   tau = ROOT.TLorentzVector()
   tau.SetPtEtaPhiM(eclean[0].pt, eclean[0].eta, eclean[0].phi, eclean[0].mass)

   if len(s_j) > 1 :
       j0 = ROOT.TLorentzVector()
       j0.SetPtEtaPhiM(s_j[0].pt, s_j[0].eta, s_j[0].phi, s_j[0].mass)
       if tau.DeltaR(j0) < 0.2 :
           j = ROOT.TLorentzVector()
           j.SetPtEtaPhiM(s_j[1].pt, s_j[1].eta, s_j[1].phi, s_j[1].mass)
       else:
           j = j0
   if len(s_j) == 1:
       j = ROOT.TLorentzVector()
       j.SetPtEtaPhiM(s_j[0].pt, s_j[0].eta, s_j[0].phi, s_j[0].mass)

   m = ROOT.TLorentzVector()
   m.SetPtEtaPhiM(met_pt, 0, met_phi, 0)

   if (e+tau).M() > 1:

       h['hETauBaseline'].Fill((e+tau).M(), genweight)

       h['hETau_JetPt'].Fill(j.Pt(), genweight)
       h['hETau_MetPt'].Fill(met_pt, genweight)
       h['hETau_ElectronPt'].Fill(met_pt, genweight)

       if e.DeltaR(tau) < 0.4 and e.DeltaR(j) > 0.8 and tau.DeltaR(j) > 0.8 :
           if met_pt > 100:
               if abs(m.DeltaPhi(e)) < 1 and abs(m.DeltaPhi(j)) > 2:
                   if len(s_isomu) == 0:
                       h['hETauSelection_veto'].Fill((e+tau).M(), genweight)

                   if e.DeltaR(tau) > 0.05 :
                       h['hETauSelection'].Fill((e+tau).M(), genweight)
                       h['hETauSelection_JetPt'].Fill(j.Pt(), genweight)
                       h['hETauSelection_MetPt'].Fill(met_pt, genweight)
                       h['hETauSelection_ElectronPt'].Fill(e.Pt(), genweight)

                   if ( j.Pt() > 500 and ( isSingleJet == 1 or isHT == 1 ) ) \
                      or ( e.Pt() > 35 and (isIsoEle == 1 ) ) \
                      or ( isHTMHT == 1 ) :
                       h['hETauTrigger_JetPt_MetPt'].Fill(j.Pt(), met_pt, genweight)


                   if ( j.Pt() > 500 and ( isSingleJet == 1 or isHT == 1 ) ) \
                      or ( e.Pt() > 35 and (isIsoEle == 1 ) ) \
                      or ( met_pt > 200 and isHTMHT == 1 ) :
                       if e.DeltaR(tau) > 0.05 :
                           h['hETauTrigger_JetPt'].Fill(j.Pt(), genweight)
                           h['hETauTrigger_MetPt'].Fill(met_pt, genweight)

                   if ( j.Pt() > 500 and ( isSingleJet == 1 or isHT == 1 ) ) \
                      or ( e.Pt() > 35 and (isIsoEle == 1 ) ) \
                      or ( j.Pt() > 100 and met_pt > 200 and isHTMHT == 1 ) :
                       h['hETauTrigger_MetHT_Loose'].Fill((e+tau).M(), genweight)
                       if e.DeltaR(tau) > 0.05 :
                           h['hETauTrigger_MetHT_Loose_wdR'].Fill((e+tau).M(), genweight)

              #         isETauLoose = 1

                   if ( j.Pt() > 500 and ( isSingleJet == 1 or isHT == 1 ) ) \
                      or ( e.Pt() > 35 and (isIsoEle == 1 ) ) \
                      or ( j.Pt() > 500 and met_pt > 200 and isHTMHT == 1 ) :
                       h['hETauTrigger_MetHT_Tight'].Fill((e+tau).M(), genweight)
#                       isETauTight = 1
                       h['hTauEvents'].Fill(5,1)
                       h['hETauTrigger_dRl'].Fill((e+tau).M(),e.DeltaR(tau), genweight)
                       if e.DeltaR(tau) > 0.05 :
                           isETauTight = 1
                           h['hETauTrigger_MetHT_Tight_wdR'].Fill((e+tau).M(), genweight)


                       if len(s_e) > 1:
                           e2 = ROOT.TLorentzVector()
                           e2.SetPtEtaPhiM(s_e[1].pt, s_e[1].eta, s_e[1].phi, s_e[1].mass)
                           if e.DeltaR(e2) < 0.4 and e2.DeltaR(j) > 0.8:
                               h['hETauTrigger_EE'].Fill((e+e2).M(), genweight)
                               h['hETauTrigger_EE_ETau'].Fill((e+tau).M(), genweight)

   return isETauTight

def TauTau_Channel():

   isTauTau = 0

   tau1 = ROOT.TLorentzVector()
   tau1.SetPtEtaPhiM(boosted[0].pt, boosted[0].eta, boosted[0].phi, boosted[0].mass)

   tau2 = ROOT.TLorentzVector()
   tau2.SetPtEtaPhiM(boosted[1].pt, boosted[1].eta, boosted[1].phi, boosted[1].mass)

   if len(s_j) > 1:
       j0 = ROOT.TLorentzVector()
       j0.SetPtEtaPhiM(s_j[0].pt, s_j[0].eta, s_j[0].phi, s_j[0].mass)
       if ( j0.DeltaR(tau1) < 0.2 ) or ( j0.DeltaR(tau2) < 0.2 ) :
           j = ROOT.TLorentzVector()
           j.SetPtEtaPhiM(s_j[1].pt, s_j[1].eta, s_j[1].phi, s_j[1].mass)
       else :
           j = j0

   if len(s_j) == 1 :
       j = ROOT.TLorentzVector()
       j.SetPtEtaPhiM(s_j[0].pt, s_j[0].eta, s_j[0].phi, s_j[0].mass)

   m = ROOT.TLorentzVector()
   m.SetPtEtaPhiM(met_pt, 0, met_phi, 0)

   if (tau1+tau2).M() > 1:

       h['hTauTauBaseline'].Fill((tau1+tau2).M(), genweight)

       h['hTauTau_JetPt'].Fill(j.Pt(), genweight)
       h['hTauTau_MetPt'].Fill(met_pt, genweight)


       if tau1.DeltaR(tau2) < 0.4 and tau1.DeltaR(j) > 0.8 and tau2.DeltaR(j) > 0.8:
           if met_pt > 100:
               if len(s_isomu) == 0 and len(s_isoe) == 0:
                   h['hTauTauSelection_veto'].Fill((tau1+tau2).M(), genweight)

               h['hTauTauSelection'].Fill((tau1+tau2).M(), genweight)
               h['hTauTauSelection_JetPt'].Fill(j.Pt(), genweight)
               h['hTauTauSelection_MetPt'].Fill(met_pt, genweight)

               if ( j.Pt() > 500 and ( isSingleJet == 1 or isHT == 1 ) ) \
                  or ( isHTMHT == 1 ) :
                   h['hTauTauTrigger_JetPt_MetPt'].Fill(j.Pt(), met_pt, genweight)

               if ( j.Pt() > 500 and ( isSingleJet == 1 or isHT == 1 ) ) \
                  or ( met_pt > 200 and isHTMHT == 1 ) :
                   h['hTauTauTrigger_JetPt'].Fill(j.Pt(), genweight)
                   h['hTauTauTrigger_MetPt'].Fill(met_pt, genweight)

               if ( j.Pt() > 500 and ( isSingleJet == 1 or isHT == 1 ) ) \
                  or ( j.Pt() > 100 and met_pt > 200 and isHTMHT == 1 ) :
#                   if isETauLoose == 0:
                   h['hTauTauTrigger_MetHT_Loose'].Fill((tau1+tau2).M(), genweight)

               if ( j.Pt() > 500 and ( isSingleJet == 1 or isHT == 1 ) ) \
                  or ( j.Pt() > 500 and met_pt > 200 and isHTMHT == 1 ) :
                   #if isETauTight == 0: 
                   h['hTauTauTrigger_MetHT_Tight'].Fill((tau1+tau2).M(), genweight)
                   h['hTauEvents'].Fill(6,1)

                   h['hTauTauTrigger_MetHT_nLepton'].Fill(0.5, len(s_e))
                   h['hTauTauTrigger_MetHT_nLepton'].Fill(1.5, len(s_mu))

                   if len(s_e) > 1 :
                       e = ROOT.TLorentzVector()
                       e.SetPtEtaPhiM(s_e[0].pt, s_e[0].eta, s_e[0].phi, s_e[0].mass)
                       e2 = ROOT.TLorentzVector()
                       e2.SetPtEtaPhiM(s_e[1].pt, s_e[1].eta, s_e[1].phi, s_e[1].mass)

                       if e.DeltaR(e2) < 0.4 and e.DeltaR(j) > 0.8 and e2.DeltaR(j) > 0.8 :
                           h['hTauTauTrigger_EE'].Fill((e+e2).M(), genweight)
                           h['hTauTauTrigger_EE_TauTau'].Fill((tau1+tau2).M(), genweight)

                   if len(s_mu) > 1 :
                       mu = ROOT.TLorentzVector()
                       mu.SetPtEtaPhiM(s_mu[0].pt, s_mu[0].eta, s_mu[0].phi, s_mu[0].mass)
                       mu2 = ROOT.TLorentzVector()
                       mu2.SetPtEtaPhiM(s_mu[1].pt, s_mu[1].eta, s_mu[1].phi, s_mu[1].mass)

                       if mu.DeltaR(mu2) < 0.4 and mu.DeltaR(j) > 0.8 and mu2.DeltaR(j) > 0.8 :
                           h['hTauTauTrigger_MuMu'].Fill((mu+mu2).M(), genweight)
                           h['hTauTauTrigger_MuMu_TauTau'].Fill((tau1+tau2).M(), genweight)

   return isTauTau


inputFileNames=open(inputFileList, 'r')
for inputFileName in inputFileNames:
    inputFileName=inputFileName.replace("\n","")
    print(inputFileName.replace("\n",""))

    fchain.Add(inputFileName)
    fchain.AddFriend('tcpTrigNtuples/triggerTree', inputFileName)
    fchain.AddFriend('lumiSummary/lumiTree', inputFileName)

    jets = ROOT.JetInfoDS()
    muons = ROOT.MuonInfoDS()
    electrons = ROOT.ElectronInfoDS()
    tausUnCleaned = ROOT.TauInfoDS()
    tausECleaned = ROOT.TauInfoDS()
    tausMCleaned = ROOT.TauInfoDS()
#    tausLowPtECleaned = ROOT.TauInfoDS()
    tausBoosted = ROOT.TauInfoDS()

    fchain.SetBranchAddress("Jets", ROOT.AddressOf(jets))
    fchain.SetBranchAddress("Muons", ROOT.AddressOf(muons))
    fchain.SetBranchAddress("Electrons", ROOT.AddressOf(electrons))
    fchain.SetBranchAddress("TausUnCleaned", ROOT.AddressOf(tausUnCleaned))
    fchain.SetBranchAddress("TausECleaned", ROOT.AddressOf(tausECleaned))
    fchain.SetBranchAddress("TausMCleaned", ROOT.AddressOf(tausMCleaned))
#    fchain.SetBranchAddress("TausLowPtECleaned", ROOT.AddressOf(tausLowPtECleaned))
    fchain.SetBranchAddress("TausBoosted", ROOT.AddressOf(tausBoosted))

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
   s_isomu = []
   s_isoe = []

   unclean = []
   eclean = []
   mclean = []
#       lowclean = []
   boosted = []

   if tausUnCleaned.size()>0:
       for i in range(tausUnCleaned.size()):
           tau = tausUnCleaned.at(i)
           if tau.mvaid >= 4:
               h['hTauUnCleanedPt'].Fill(tau.pt, genweight)
               unclean+=[tau]

   if tausECleaned.size()>0:
       for i in range(tausECleaned.size()):
           tau = tausECleaned.at(i)
           if tau.mvaid >= 4:
               h['hECleanedPt'].Fill(tau.pt, genweight)
               eclean+=[tau]

   if tausMCleaned.size()>0:
       for i in range(tausMCleaned.size()):
           tau = tausMCleaned.at(i)
           if tau.mvaid >= 4:
               h['hMuCleanedPt'].Fill(tau.pt, genweight)
               mclean+=[tau]

   # if tausLowPtECleaned.size()>0:
   #     for i in range(tausLowPtECleaned.size()):
   #         tau = tausLowPtECleaned.at(i)
   #         if tau.mvaid >= 4:
   #             lowclean+=[tau]

   if tausBoosted.size()>0:
       for i in range(tausBoosted.size()):
           tau = tausBoosted.at(i)
           if tau.mvaid >= 1:
               h['hTauBoostedPt'].Fill(tau.pt, genweight)
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
             if muon.iso < 0.25:
                 h['hIsoMuonPt'].Fill(muon.pt)
                 s_isomu+=[muon]

   if electrons.size()>0:
      for i in range(electrons.size()):
         electron = electrons.at(i)
         if electron.id >= 1 :
             h['hElectronPt'].Fill(electron.pt)
             s_e+=[electron]
             if electron.iso >= 1:
                 h['hIsoElectronPt'].Fill(electron.pt)
                 s_isoe+=[electron]

   s_j.sort(key=lambda x: x.pt, reverse=True)
   s_e.sort(key=lambda x: x.pt, reverse=True)
   s_isoe.sort(key=lambda x: x.pt, reverse=True)
   s_mu.sort(key=lambda x: x.pt, reverse=True)
   s_isomu.sort(key=lambda x: x.pt, reverse=True)
   unclean.sort(key=lambda x: x.pt, reverse=True)
   eclean.sort(key=lambda x: x.pt, reverse=True)
   mclean.sort(key=lambda x: x.pt, reverse=True)
#       lowclean.sort(key=lambda x: x.pt, reverse=True)
   boosted.sort(key=lambda x: x.pt, reverse=True)


   if len(s_mu) > 1 and len(s_j) > 0 and len(s_b) == 0 and s_mu[0].charge*s_mu[1].charge < 0 : 
       if MuMu_Channel() == 1: continue

   if len(s_mu) > 0 and len(s_e) > 0 and len(s_j) > 0 and len(s_b) == 0 and s_e[0].charge*s_mu[0].charge < 0 : 
       if EMu_Channel() == 1: continue

   if len(s_e) > 1 and len(s_j) > 0 and len(s_b) == 0 and s_e[0].charge*s_e[1].charge < 0 : 
       if EE_Channel() == 1: continue

   if len(s_mu) > 0 and len(mclean) > 0 and len(s_j) > 0 and len(s_b) == 0 and mclean[0].charge*s_mu[0].charge < 0 : 
       if MuTau_Channel() == 1 : continue

   if len(s_e) > 0 and len(eclean) > 0 and len(s_j) > 0 and len(s_b) == 0 and eclean[0].charge*s_e[0].charge < 0 : 
       if ETau_Channel() == 1 : continue

   if len(boosted) > 1 and len(s_j) > 0 and len(s_b) == 0 and boosted[0].charge*boosted[1].charge < 0 : 
       TauTau_Channel()


out.cd()

for key in h.keys():
    h[key].Write()

out.Close()

print("--- %s seconds ---" % (time.time() - start_time))
