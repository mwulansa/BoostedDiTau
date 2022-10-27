import ROOT, sys, os
import numpy as np
import time

start_time = time.time()

opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]

outputTitle = "h_studyEMuDataMC"

isNonIso = False
triggerStudy = False
bEffStudy = False

isMuonEGDataset = 0
isSingleMuonDataset = 0

if "-b" in opts:
    isData = 0

if "-de" in opts :
    isData = 1
    isMuonEGDataset = 1
    isSingleMuonDataset = 0

if "-ds" in opts :
    isData = 1
    isMuonEGDataset = 0
    isSingleMuonDataset = 1

ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/JetInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/MuonInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/ElectronInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/TauInfoDS.h"')

inputFileListName=sys.argv[1]
inputFileList=inputFileListName

if len(sys.argv)>2 and not sys.argv[2].startswith("-"):
    outputFileDir=sys.argv[2]
else:
    outputFileDir = "./output/"

if isNonIso == False:
    outputFileName = outputFileDir+outputTitle+"_Isolated_"+inputFileListName.split("/")[-1].replace(".txt",".root")
if isNonIso == True:
    outputFileName = outputFileDir+outputTitle+"_Non-Isolated_"+inputFileListName.split("/")[-1].replace(".txt",".root")

out=ROOT.TFile.Open(outputFileName,'recreate')
print(outputFileName)

fchain = ROOT.TChain('tcpNtuples/analysisTree')
chain2 = ROOT.TChain('tcpTrigNtuples/triggerTree')
if isData == 0:
    chain3 = ROOT.TChain('lumiSummary/lumiTree')
    chain4 = ROOT.TChain('tcpGenNtuples/genTree')

pi = np.pi

h = {}

event_cut = {
    'jetPt': 70,
    'dRl': 0.4,
    'dRltau': 0.05,
    'dRlj': 0.8,
    'metcut': 100.0,
    'mtcut': 50.0,
    'dPhiml': 1,
    'dPhimj': 2,
    'mass' : 5
}


def define_event_histogram(region):

    h[region+"_Count"] = ROOT.TH1F (region+"_Count", region+"_Count ; Events ; Events ", 1, 0, 1)

    h[region+"_Mass"] = ROOT.TH1F (region+"_Mass", region+"_Mass ; M_{vis.} (GeV) ; Events ", 100, 0, 100)
    h[region+"_Lepton1Pt"] = ROOT.TH1F (region+"_Lepton1Pt", region+"_Lepton1Pt ; P_{T} (GeV) ; Events ", 500, 0, 500)
    h[region+"_Lepton2Pt"] = ROOT.TH1F (region+"_Lepton2Pt", region+"_Lepton2Pt ; P_{T} (GeV) ; Events ", 500, 0, 500)
    h[region+"_JetPt"] = ROOT.TH1F (region+"_JetPt", region+"_JetPt ; JetP_{T} (GeV) ; Events ", 2000, 0, 2000)
    h[region+"_MetPt"] = ROOT.TH1F (region+"_MetPt", region+"_MetPt ; MET (GeV) ; Events ", 500, 0, 500)
    h[region+"_Nj"] = ROOT.TH1F (region+"_Nj", region+"_Nj ; N_{j} ; Events ", 10, 0, 10)
    h[region+"_dRl"] = ROOT.TH1F (region+"_dRl", region+"_dRl ; dR(leptons) ; Events", 100, 0, 5)
    h[region+"_dRj"] = ROOT.TH1F (region+"_dRj", region+"_dRj ; dR(jet, ditau) ; Events", 100, 0, 5)
    h[region+"_dPhil"] = ROOT.TH2F (region+"_dPhil", region+"_dPhil ; dPhi(met,lepton1) ; dPhi(met,lepton2)",  100, -pi, pi, 100, -pi, pi)
    h[region+"_dPhi"] = ROOT.TH2F (region+"_dPhi", region+"_dPhi ; dPhi(met,ditau) ; dPhi(met,jet)",  100, -pi, pi, 100, -pi, pi)

    h[region+"_Mtl1"] = ROOT.TH1F (region+"_Mtl1", region+"_Mtl1 ; m_{T}(met, lepton1) (GeV) ; Events ", 150, 0, 150)
    h[region+"_Mtl2"] = ROOT.TH1F (region+"_Mtl2", region+"_Mtl2 ; m_{T}(met, lepton2) (GeV) ; Events ", 150, 0, 150)
    h[region+"_Mtl"] = ROOT.TH1F (region+"_Mtl", region+"_Mtl ; m_{T}(met, leptons) (GeV) ; Events ", 150, 0, 150)

    h[region+"_cosMl1"] = ROOT.TH1F (region+"_cosMl1", region+"_cosMl1 ; cos(met, l1) ; Events ", 100, -1, 1)
    h[region+"_cosMl2"] = ROOT.TH1F (region+"_cosMl2", region+"_cosMl2 ; cos(met, l2) ; Events ", 100, -1, 1)
    h[region+"_cosl"] = ROOT.TH1F (region+"_cosl", region+"_cosl ; cos(l1, l2) ; Events ", 100, -1, 1)

    h[region+"_bJetPt"] = ROOT.TH1F (region+"_bJetPt", region+"_bJetPt ; P_{T} (GeV) ; Events ", 2000, 0, 2000)

    if region in isoPlotRegions:
        h[region+"_muIso"] = ROOT.TH1F (region+"_muIso", region+"_muIso ; muIso ; Events ", 100, 0, 1)


def define_b_eff_histogram(region):

    h[region+"_BFlavour_JetPt"] = ROOT.TH1F (region+"_BFlavour_JetPt", region+"_BFlavour_JetPt ; P_{T} (GeV) ; Events", 2000, 0 , 2000)
    h[region+"_BFlavour_JetEta"] = ROOT.TH1F (region+"_BFlavour_JetEta", region+"_BFlavour_JetEta ; Eta ; Events", 100, 0, 2.5)
    
    h[region+"_BFlavour_JetPtEta"] = ROOT.TH2F (region+"_BFlavour_JetPtEta", region+"_BFlavour_JetPtEta ; P_{T} (GeV) ; Eta", 2000, 0, 2000, 100, 0, 2.5)
    
    h[region+"_BFlavour_BTagged_JetPt"] = ROOT.TH1F (region+"_BFlavour_BTagged_JetPt", region+"_BFlavour_BTagged_JetPt ; P_{T} (GeV) ; Events", 2000, 0 , 2000)
    h[region+"_BFlavour_BTagged_JetEta"] = ROOT.TH1F (region+"_BFlavour_BTagged_JetEta", region+"_BFlavour_BTagged_JetEta ; Eta ; Events", 100, 0, 2.5)

    h[region+"_BFlavour_BTagged_JetPtEta"] = ROOT.TH2F (region+"_BFlavour_BTagged_JetPtEta", region+"_BFlavour_BTagged_JetPtEta ; P_{T} (GeV) ; Eta", 2000, 0, 2000, 100, 0, 2.5)

    h[region+"_CFlavour_JetPt"] = ROOT.TH1F (region+"_CFlavour_JetPt", region+"_CFlavour_JetPt ; P_{T} (GeV) ; Events", 2000, 0 , 2000)
    h[region+"_CFlavour_JetEta"] = ROOT.TH1F (region+"_CFlavour_JetEta", region+"_CFlavour_JetEta ; Eta ; Events", 100, 0 , 2.5)

    h[region+"_CFlavour_JetPtEta"] = ROOT.TH2F (region+"_CFlavour_JetPtEta", region+"_CFlavour_JetPtEta ; P_{T} (GeV) ; Eta", 2000, 0, 2000, 100, 0, 2.5)

    h[region+"_CFlavour_BTagged_JetPt"] = ROOT.TH1F (region+"_CFlavour_BTagged_JetPt", region+"_CFlavour_BTagged_JetPt ; P_{T} (GeV) ; Events", 2000, 0 , 2000)
    h[region+"_CFlavour_BTagged_JetEta"] = ROOT.TH1F (region+"_CFlavour_BTagged_JetEta", region+"_CFlavour_BTagged_JetEta ; Eta ; Events", 100, 0 , 2.5)

    h[region+"_CFlavour_BTagged_JetPtEta"] = ROOT.TH2F (region+"_CFlavour_BTagged_JetPtEta", region+"_CFlavour_BTagged_JetPtEta ; P_{T} (GeV) ; Eta", 2000, 0, 2000, 100, 0, 2.5)

    h[region+"_LFlavour_JetPt"] = ROOT.TH1F (region+"_LFlavour_JetPt", region+"_LFlavour_JetPt ; P_{T} (GeV) ; Events", 2000, 0 , 2000)
    h[region+"_LFlavour_JetEta"] = ROOT.TH1F (region+"_LFlavour_JetEta", region+"_LFlavour_JetEta ; Eta ; Events", 100, 0 , 2.5)

    h[region+"_LFlavour_JetPtEta"] = ROOT.TH2F (region+"_LFlavour_JetPtEta", region+"_LFlavour_JetPtEta ; P_{T} (GeV) ; Eta", 2000, 0, 2000, 100, 0, 2.5)

    h[region+"_LFlavour_BTagged_JetPt"] = ROOT.TH1F (region+"_LFlavour_BTagged_JetPt", region+"_LFlavour_BTagged_JetPt ; P_{T} (GeV) ; Events", 2000, 0 , 2000)
    h[region+"_LFlavour_BTagged_JetEta"] = ROOT.TH1F (region+"_LFlavour_BTagged_JetEta", region+"_LFlavour_BTagged_JetEta ; Eta ; Events", 100, 0 , 2.5)

    h[region+"_LFlavour_BTagged_JetPtEta"] = ROOT.TH2F (region+"_LFlavour_BTagged_JetPtEta", region+"_LFlavour_BTagged_JetPtEta ; P_{T} (GeV) ; Eta", 2000, 0, 2000, 100, 0, 2.5)

    h[region+"_Loosebjet_Pt"] = ROOT.TH1F (region+"_Loosebjet_Pt", region+"_Loosebjet_Pt ; Pt (GeV) ; Events", 2000, 0, 2000)
    h[region+"_Loosebjet_Eta"] = ROOT.TH1F (region+"_Loosebjet_Eta", region+"_Loosebjet_Eta ; Eta ; Events", 100, 0, 2.5)

    h[region+"_Mediumbjet_Pt"] = ROOT.TH1F (region+"_Mediumbjet_Pt", region+"_Mediumbjet_Pt ; Pt (GeV) ; Events", 2000, 0, 2000)
    h[region+"_Mediumbjet_Eta"] = ROOT.TH1F (region+"_Mediumbjet_Eta", region+"_Mediumbjet_Eta ; Eta ; Events", 100, 0, 2.5)

    h[region+"_Tightbjet_Pt"] = ROOT.TH1F (region+"_Tightbjet_Pt", region+"_Tightbjet_Pt ; Pt (GeV) ; Events", 2000, 0, 2000)
    h[region+"_Tightbjet_Eta"] = ROOT.TH1F (region+"_Tightbjet_Eta", region+"_Tightbjet_Eta ; Eta ; Events", 100, 0, 2.5)


def define_bjet_histogram(region):

    h[region+'_dRbjet'] = ROOT.TH1F (region+"_dRbjet", region+"_dRbjet ; dR(bjets) ; Events", 100, 0, 5)
    h[region+'_dRebjet1'] = ROOT.TH1F (region+"_dRebjet1", region+"_dRebjet1 ; dR(e,bjet1) ; Events", 100, 0, 5)
    h[region+'_dRebjet2'] = ROOT.TH1F (region+"_dRebjet2", region+"_dRebjet2 ; dR(e,bjet2) ; Events", 100, 0, 5)
    h[region+'_dRmubjet1'] = ROOT.TH1F (region+"_dRmubjet1", region+"_dRmubjet1 ; dR(mu,bjet1) ; Events", 100, 0, 5)
    h[region+'_dRmubjet2'] = ROOT.TH1F (region+"_dRmubjet2", region+"_dRmubjet2 ; dR(mu,bjet2) ; Events", 100, 0, 5)
    h[region+'_deepjet'] = ROOT.TH1F (region+"_deepjet", region+"_deepjet ; deepjet score ; Events", 100, 0, 1)


def define_general_histogram():

#-----Trigger bits-----

    h['isSingleJet'] = ROOT.TH1F("isSingleJet", "isSingleJet ; isSingleJet ; N", 4,-1.5,2.5)
    h['isHT'] = ROOT.TH1F("isHT", "isHT ; isHT ; N", 4,-1.5,2.5)

#-----Event counts-----

    h['hEvents'] = ROOT.TH1F ("NEvents", "Number of Events; ;N", 2, 0, 2)

    h['hMuMu_Events'] = ROOT.TH1F ("hMuMu_Events", "hMuMu_Events;;N", 6, 1, 7)
    h['hEMu_Events'] = ROOT.TH1F ("hEMu_Events", "hEMu_Events;;N", 6, 1, 7)
    h['hEMu_nonIsoMu_Events'] = ROOT.TH1F ("hEMu_nonIsoMu_Events", "hEMu_Events;;N", 6, 1, 7)
    h['hEMu_nonIsoE_Events'] = ROOT.TH1F ("hEMu_nonIsoE_Events", "hEMu_Events;;N", 6, 1, 7)
    h['hEMu_nonIso_Events'] = ROOT.TH1F ("hEMu_nonIso_Events", "hEMu_Events;;N", 6, 1, 7)
    h['hEMu_withoutMuMu_Events'] = ROOT.TH1F ("hEMu_withoutMuMu_Events", "hEMu_Events;;N", 6, 1, 7)
    h['hEMu_checkMuMu_Events'] = ROOT.TH1F ("hEMu_checkMuMu_Events", "hEMu_Events;;N", 6, 1, 7)

#-----Bjets-----

    h["kin_baseline_deepjet"] = ROOT.TH1F ("kin_baseline_deepjet", "kin_baseline_deepjet ; deepjet score ; Events ", 100, 0, 1)
    h["kin_ePtcut_deepjet"] = ROOT.TH1F ("kin_ePtcut_deepjet", "kin_ePtcut_deepjet ; deepjet score ; Events ", 100, 0, 1)
    h["kin_lPtcut_deepjet"] = ROOT.TH1F ("kin_lPtcut_deepjet", "kin_lPtcut_deepjet ; deepjet score ; Events ", 100, 0, 1)
    h["kin_jetPtcut_deepjet"] = ROOT.TH1F ("kin_jetPtcut_deepjet", "kin_jetPtcut_deepjet ; deepjet score ; Events ", 100, 0, 1)

    h["hEMu_Nbjet"] = ROOT.TH1F ("hEMu_Nbjet", "hEMu_Nbjet ; N_{bjet} ; Events ", 10, 0, 10)
    h["hEMu_dRcut_Nbjet"] = ROOT.TH1F ("hEMu_dRcut_Nbjet", "hEMu_dRcut_Nbjet ; N_{bjet} ; Events ", 10, 0, 10)
    h["hEMu_lowMET_Nbjet"] = ROOT.TH1F ("hEMu_lowMET_Nbjet", "hEMu_lowMET_Nbjet ; N_{bjet} ; Events ", 10, 0, 10)
    h["hEMu_highMET_Nbjet"] = ROOT.TH1F ("hEMu_highMET_Nbjet", "hEMu_highMET_Nbjet ; N_{bjet} ; Events ", 10, 0, 10)
    h["hEMu_dRcut_lowMET_dPhicut_Nbjet"] = ROOT.TH1F ("hEMu_dRcut_lowMET_dPhicut_Nbjet", "hEMu_dRcut_lowMET_dPhicut_Nbjet ; N_{bjet} ; Events ", 10, 0, 10)
    h["hEMu_dRcut_lowMET_Nbjet"] = ROOT.TH1F ("hEMu_dRcut_lowMET_Nbjet", "hEMu_dRcut_lowMET_Nbjet ; N_{bjet} ; Events ", 10, 0, 10)
    h["hEMu_dRcut_highMET_Nbjet"] = ROOT.TH1F ("hEMu_dRcut_highMET_Nbjet", "hEMu_dRcut_highMET_Nbjet ; N_{bjet} ; Events ", 10, 0, 10)

    h["hEMu_deepjet"] = ROOT.TH1F ("hEMu_deepjet", "hEMu_deepjet ; deepjet score ; Events ", 100, 0, 1)
    h["hEMu_dRcut_deepjet"] = ROOT.TH1F ("hEMu_dRcut_deepjet", "hEMu_dRcut_deepjet ; deepjet score ; Events ", 100, 0, 1)
    h["hEMu_loosedRlcut_deepjet"] = ROOT.TH1F ("hEMu_loosedRlcut_deepjet", "hEMu_loosedRlcut_deepjet ; deepjet score ; Events ", 100, 0, 1)
    h["hEMu_nodRlcut_deepjet"] = ROOT.TH1F ("hEMu_nodRlcut_deepjet", "hEMu_nodRlcut_deepjet ; deepjet score ; Events ", 100, 0, 1)
    h["hEMu_lowMET_deepjet"] = ROOT.TH1F ("hEMu_lowMET_deepjet", "hEMu_lowMET_deepjet ; deepjet score ; Events ", 100, 0, 1)
    h["hEMu_highMET_deepjet"] = ROOT.TH1F ("hEMu_highMET_deepjet", "hEMu_highMET_deepjet ; deepjet score ; Events ", 100, 0, 1)
    h["hEMu_dRcut_lowMET_deepjet"] = ROOT.TH1F ("hEMu_dRcut_lowMET_deepjet", "hEMu_dRcut_lowMET_deepjet ; deepjet score ; Events ", 100, 0, 1)
    h["hEMu_dRcut_highMET_deepjet"] = ROOT.TH1F ("hEMu_dRcut_highMET_deepjet", "hEMu_dRcut_lowMET_deepjet ; deepjet score ; Events ", 100, 0, 1)
    h["hEMu_dRcut_lowMET_dPhicut_deepjet"] = ROOT.TH1F ("hEMu_dRcut_lowMET_dPhicut_deepjet", "hEMu_dRcut_lowMET_dPhicut_deepjet ; deepjet score ; Events ", 100, 0, 1)

    h["hEMu_checkMuMu_Nbjet"] = ROOT.TH1F ("hEMu_checkMuMu_Nbjet", "hEMu_checkMuMu_Nbjet ; N_{bjet} ; Events ", 10, 0, 10)
    h["hEMu_checkMuMu_dRcut_Nbjet"] = ROOT.TH1F ("hEMu_checkMuMu_dRcut_Nbjet", "hEMu_checkMuMu_dRcut_Nbjet ; N_{bjet} ; Events ", 10, 0, 10)
    h["hEMu_checkMuMu_lowMET_Nbjet"] = ROOT.TH1F ("hEMu_checkMuMu_lowMET_Nbjet", "hEMu_checkMuMu_lowMET_Nbjet ; N_{bjet} ; Events ", 10, 0, 10)
    h["hEMu_checkMuMu_highMET_Nbjet"] = ROOT.TH1F ("hEMu_checkMuMu_highMET_Nbjet", "hEMu_checkMuMu_highMET_Nbjet ; N_{bjet} ; Events ", 10, 0, 10)
    h["hEMu_checkMuMu_dRcut_lowMET_dPhicut_Nbjet"] = ROOT.TH1F ("hEMu_checkMuMu_dRcut_lowMET_dPhicut_Nbjet", "hEMu_checkMuMu_dRcut_lowMET_dPhicut_Nbjet ; N_{bjet} ; Events ", 10, 0, 10)
    h["hEMu_checkMuMu_dRcut_lowMET_Nbjet"] = ROOT.TH1F ("hEMu_checkMuMu_dRcut_lowMET_Nbjet", "hEMu_checkMuMu_dRcut_lowMET_Nbjet ; N_{bjet} ; Events ", 10, 0, 10)
    h["hEMu_checkMuMu_dRcut_highMET_Nbjet"] = ROOT.TH1F ("hEMu_checkMuMu_dRcut_highMET_Nbjet", "hEMu_checkMuMu_dRcut_highMET_Nbjet ; N_{bjet} ; Events ", 10, 0, 10)

    h["hEMu_checkMuMu_deepjet"] = ROOT.TH1F ("hEMu_checkMuMu_deepjet", "hEMu_checkMuMu_deepjet ; deepjet score ; Events ", 100, 0, 1)
    h["hEMu_checkMuMu_dRcut_deepjet"] = ROOT.TH1F ("hEMu_checkMuMu_dRcut_deepjet", "hEMu_checkMuMu_dRcut_deepjet ; deepjet score ; Events ", 100, 0, 1)
    h["hEMu_checkMuMu_loosedRlcut_deepjet"] = ROOT.TH1F ("hEMu_checkMuMu_loosedRlcut_deepjet", "hEMu_checkMuMu_loosedRlcut_deepjet ; deepjet score ; Events ", 100, 0, 1)
    h["hEMu_checkMuMu_nodRlcut_deepjet"] = ROOT.TH1F ("hEMu_checkMuMu_nodRlcut_deepjet", "hEMu_checkMuMu_nodRlcut_deepjet ; deepjet score ; Events ", 100, 0, 1)
    h["hEMu_checkMuMu_lowMET_deepjet"] = ROOT.TH1F ("hEMu_checkMuMu_lowMET_deepjet", "hEMu_checkMuMu_lowMET_deepjet ; deepjet score ; Events ", 100, 0, 1)
    h["hEMu_checkMuMu_highMET_deepjet"] = ROOT.TH1F ("hEMu_checkMuMu_highMET_deepjet", "hEMu_checkMuMu_highMET_deepjet ; deepjet score ; Events ", 100, 0, 1)
    h["hEMu_checkMuMu_dRcut_lowMET_deepjet"] = ROOT.TH1F ("hEMu_checkMuMu_dRcut_lowMET_deepjet", "hEMu_checkMuMu_dRcut_lowMET_deepjet ; deepjet score ; Events ", 100, 0, 1)
    h["hEMu_checkMuMu_dRcut_highMET_deepjet"] = ROOT.TH1F ("hEMu_checkMuMu_dRcut_highMET_deepjet", "hEMu_checkMuMu_dRcut_lowMET_deepjet ; deepjet score ; Events ", 100, 0, 1)
    h["hEMu_checkMuMu_dRcut_lowMET_dPhicut_deepjet"] = ROOT.TH1F ("hEMu_checkMuMu_dRcut_lowMET_dPhicut_deepjet", "hEMu_checkMuMu_dRcut_lowMET_dPhicut_deepjet ; deepjet score ; Events ", 100, 0, 1)

    h["hMuMu_Nbjet"] = ROOT.TH1F ("hMuMu_Nbjet", "hMuMu_Nbjet ; N_{bjet} ; Events ", 10, 0, 10)
    h["hMuMu_dRcut_Nbjet"] = ROOT.TH1F ("hMuMu_dRcut_Nbjet", "hMuMu_dRcut_Nbjet ; N_{bjet} ; Events ", 10, 0, 10)
    h["hMuMu_lowMET_Nbjet"] = ROOT.TH1F ("hMuMu_lowMET_Nbjet", "hMuMu_lowMET_Nbjet ; N_{bjet} ; Events ", 10, 0, 10)
    h["hMuMu_highMET_Nbjet"] = ROOT.TH1F ("hMuMu_highMET_Nbjet", "hMuMu_highMET_Nbjet ; N_{bjet} ; Events ", 10, 0, 10)
    h["hMuMu_dRcut_lowMET_dPhicut_Nbjet"] = ROOT.TH1F ("hMuMu_dRcut_lowMET_dPhicut_Nbjet", "hMuMu_dRcut_lowMET_dPhicut_Nbjet ; N_{bjet} ; Events ", 10, 0, 10)
    h["hMuMu_dRcut_lowMET_Nbjet"] = ROOT.TH1F ("hMuMu_dRcut_lowMET_Nbjet", "hMuMu_dRcut_lowMET_Nbjet ; N_{bjet} ; Events ", 10, 0, 10)
    h["hMuMu_dRcut_highMET_Nbjet"] = ROOT.TH1F ("hMuMu_dRcut_highMET_Nbjet", "hMuMu_dRcut_highMET_Nbjet ; N_{bjet} ; Events ", 10, 0, 10)

    h["hMuMu_deepjet"] = ROOT.TH1F ("hMuMu_deepjet", "hMuMu_deepjet ; deepjet score ; Events ", 100, 0, 1)
    h["hMuMu_dRcut_deepjet"] = ROOT.TH1F ("hMuMu_dRcut_deepjet", "hMuMu_dRcut_deepjet ; deepjet score ; Events ", 100, 0, 1)
    h["hMuMu_lowMET_deepjet"] = ROOT.TH1F ("hMuMu_lowMET_deepjet", "hMuMu_lowMET_deepjet ; deepjet score ; Events ", 100, 0, 1)
    h["hMuMu_highMET_deepjet"] = ROOT.TH1F ("hMuMu_highMET_deepjet", "hMuMu_highMET_deepjet ; deepjet score ; Events ", 100, 0, 1)
    h["hMuMu_dRcut_lowMET_deepjet"] = ROOT.TH1F ("hMuMu_dRcut_lowMET_deepjet", "hMuMu_dRcut_lowMET_deepjet ; deepjet score ; Events ", 100, 0, 1)
    h["hMuMu_dRcut_highMET_deepjet"] = ROOT.TH1F ("hMuMu_dRcut_highMET_deepjet", "hMuMu_dRcut_lowMET_deepjet ; deepjet score ; Events ", 100, 0, 1)
    h["hMuMu_dRcut_lowMET_dPhicut_deepjet"] = ROOT.TH1F ("hMuMu_dRcut_lowMET_dPhicut_deepjet", "hMuMu_dRcut_lowMET_dPhicut_deepjet ; deepjet score ; Events ", 100, 0, 1)

    h["hEMu_nonIsoMu_Nbjet"] = ROOT.TH1F ("hEMu_nonIsoMu_Nbjet", "hEMu_Nbjet ; N_{bjet} ; Events ", 10, 0, 10)
    h["hEMu_nonIsoE_Nbjet"] = ROOT.TH1F ("hEMu_nonIsoE_Nbjet", "hEMu_Nbjet ; N_{bjet} ; Events ", 10, 0, 10)
    h["hEMu_nonIso_Nbjet"] = ROOT.TH1F ("hEMu_nonIso_Nbjet", "hEMu_Nbjet ; N_{bjet} ; Events ", 10, 0, 10)

    h["hEMu_deepjet_dRlbj"] = ROOT.TH2F ("hEMu_deepjet_dRlbj", "hEMu_deepjet_dRlbj ; deepjet score ; dR(l,bj) ", 100, 0, 1, 100, 0, 5 )
    h["hEMu_dRcut_deepjet_dRlbj"] = ROOT.TH2F ("hEMu_dRcut_deepjet_dRlbj", "hEMu_dRcut_deepjet_dRlbj ; deepjet score ; dR(l,bj) ", 100, 0, 1, 100, 0, 5 )

    h["hEMu_lowMET_deepjet_dRlbj"] = ROOT.TH2F ("hEMu_lowMET_deepjet_dRlbj", "hEMu_lowMET_deepjet_dRlbj ; deepjet score ; dR(l,bj) ", 100, 0, 1, 100, 0, 5 )
    h["hEMu_highMET_deepjet_dRlbj"] = ROOT.TH2F ("hEMu_highMET_deepjet_dRlbj", "hEMu_highMET_deepjet_dRlbj ; deepjet score ; dR(l,bj) ", 100, 0, 1, 100, 0, 5 )
    h["hEMu_dRcut_lowMET_deepjet_dRlbj"] = ROOT.TH2F ("hEMu_dRcut_lowMET_deepjet_dRlbj", "hEMu_dRcut_lowMET_deepjet_dRlbj ; deepjet score ; dR(l,bj) ", 100, 0, 1, 100, 0, 5 )
    h["hEMu_dRcut_highMET_deepjet_dRlbj"] = ROOT.TH2F ("hEMu_dRcut_highMET_deepjet_dRlbj", "hEMu_dRcut_highMET_deepjet_dRlbj ; deepjet score ; dR(l,bj) ", 100, 0, 1, 100, 0, 5 )

    h["hEMu_checkMuMu_deepjet_dRlbj"] = ROOT.TH2F ("hEMu_checkMuMu_deepjet_dRlbj", "hEMu_checkMuMu_deepjet_dRlbj ; deepjet score ; dR(l,bj) ", 100, 0, 1, 100, 0, 5 )
    h["hEMu_checkMuMu_dRcut_deepjet_dRlbj"] = ROOT.TH2F ("hEMu_checkMuMu_dRcut_deepjet_dRlbj", "hEMu_checkMuMu_dRcut_deepjet_dRlbj ; deepjet score ; dR(l,bj) ", 100, 0, 1, 100, 0, 5 )
    h["hEMu_checkMuMu_lowMET_deepjet_dRlbj"] = ROOT.TH2F ("hEMu_checkMuMu_lowMET_deepjet_dRlbj", "hEMu_checkMuMu_lowMET_deepjet_dRlbj ; deepjet score ; dR(l,bj) ", 100, 0, 1, 100, 0, 5 )
    h["hEMu_checkMuMu_highMET_deepjet_dRlbj"] = ROOT.TH2F ("hEMu_checkMuMu_highMET_deepjet_dRlbj", "hEMu_checkMuMu_highMET_deepjet_dRlbj ; deepjet score ; dR(l,bj) ", 100, 0, 1, 100, 0, 5 )
    h["hEMu_checkMuMu_dRcut_lowMET_deepjet_dRlbj"] = ROOT.TH2F ("hEMu_checkMuMu_dRcut_lowMET_deepjet_dRlbj", "hEMu_checkMuMu_dRcut_lowMET_deepjet_dRlbj ; deepjet score ; dR(l,bj) ", 100, 0, 1, 100, 0, 5 )
    h["hEMu_checkMuMu_dRcut_highMET_deepjet_dRlbj"] = ROOT.TH2F ("hEMu_checkMuMu_dRcut_highMET_deepjet_dRlbj", "hEMu_checkMuMu_dRcut_highMET_deepjet_dRlbj ; deepjet score ; dR(l,bj) ", 100, 0, 1, 100, 0, 5 )


    h["hMuMu_deepjet_dRlbj"] = ROOT.TH2F ("hMuMu_deepjet_dRlbj", "hMuMu_deepjet_dRlbj ; deepjet score ; dR(l,bj) ", 100, 0, 1, 100, 0, 5 )
    h["hMuMu_dRcut_deepjet_dRlbj"] = ROOT.TH2F ("hMuMu_dRcut_deepjet_dRlbj", "hMuMu_dRcut_deepjet_dRlbj ; deepjet score ; dR(l,bj) ", 100, 0, 1, 100, 0, 5 )
    h["hMuMu_lowMET_deepjet_dRlbj"] = ROOT.TH2F ("hMuMu_lowMET_deepjet_dRlbj", "hMuMu_lowMET_deepjet_dRlbj ; deepjet score ; dR(l,bj) ", 100, 0, 1, 100, 0, 5 )
    h["hMuMu_highMET_deepjet_dRlbj"] = ROOT.TH2F ("hMuMu_highMET_deepjet_dRlbj", "hMuMu_highMET_deepjet_dRlbj ; deepjet score ; dR(l,bj) ", 100, 0, 1, 100, 0, 5)
    h["hMuMu_dRcut_lowMET_deepjet_dRlbj"] = ROOT.TH2F ("hMuMu_dRcut_lowMET_deepjet_dRlbj", "hMuMu_dRcut_lowMET_deepjet_dRlbj ; deepjet score ; dR(l,bj) ", 100, 0, 1, 100, 0, 5 )
    h["hMuMu_dRcut_highMET_deepjet_dRlbj"] = ROOT.TH2F ("hMuMu_dRcut_highMET_deepjet_dRlbj", "hMuMu_dRcut_highMET_deepjet_dRlbj ; deepjet score ; dR(l,bj) ", 100, 0, 1, 100, 0, 5 )


#-----Objects-----

    h['hJetPt'] = ROOT.TH1F ("hJetPt", "Jet P_{T} ; P_{T} ; N", 1000, 0, 1000)
    h['hBtag'] = ROOT.TH1F ("hBtag", "DeepJet Score; score; N", 50, 0, 1)
    h['hBJetPt'] = ROOT.TH1F ("hBJetPt", "B-tagged Jet P_{T} ; P_{T} ; N", 1000, 0, 1000)
    h['hMuonPt'] = ROOT.TH1F ("hMuPt", "Muon P_{T} ; P_{T} ; N", 500, 0, 500)
    h['hIsoMuonPt'] = ROOT.TH1F ("hIsoMuPt", "Isolated Muon P_{T} ; P_{T} ; N", 500, 0, 500)
    h['hNonIsoMuonPt'] = ROOT.TH1F ("hNonIsoMuPt", "Non-Isolated Muon P_{T} ; P_{T} ; N", 500, 0, 500)
    h['hElectronPt'] = ROOT.TH1F ("hEPt", "Electron P_{T} ; P_{T} ; N", 500, 0, 500)
    h['hIsoElectronPt'] = ROOT.TH1F ("hIsoEPt", "Isolated Electron P_{T} ; P_{T} ; N", 500, 0, 500)
    h['hNonIsoElectronPt'] = ROOT.TH1F ("hNonIsoEPt", "Non-Isolated Electron P_{T} ; P_{T} ; N", 500, 0, 500)
    h['hMetPt'] = ROOT.TH1F ("hMetPt", "MET P_{T} ; P_{T} ; N", 500, 0, 500)

    h['hTauUnCleanedPt'] = ROOT.TH1F ("hTauUnCleanedPt", "Uncleaned Tau P_{T}; P_{T}; a.u.", 500, 0, 500)
    h['hMuCleanedPt'] = ROOT.TH1F ("hMuCleanedPt", "Mu-cleaned Tau P_{T}; P_{T}; a.u.", 500, 0, 500)
    h['hMuCleanedPt_altered'] = ROOT.TH1F ("hMuCleanedPt_altered", "Mu-cleaned Tau P_{T}; P_{T}; a.u.", 500, 0, 500)
    h['hECleanedPt'] = ROOT.TH1F ("hECleanedPt", "E-cleaned Tau P_{T}; P_{T}; a.u.", 500, 0, 500)
    h['hTauBoostedPt'] = ROOT.TH1F ("hTauBoostedPt", "Boosted Tau P_{T}; P_{T}; a.u.", 500, 0, 500)


def cut_flow_hists():

    h['MuonEvent_MuonPt_SingleMuon'] = ROOT.TH1F( "MuonEvent_MuonPt_SingleMuon" , "MuonEvent_MuonPt_SingleMuon; P_{T} (GeV), Events", 500, 0, 500)
    h['MuonEvent_MetPt_SingleMuon'] = ROOT.TH1F( "MuonEvent_MetPt_SingleMuon" , "MuonEvent_MetPt_SingleMuon; P_{T} (GeV), Events", 500, 0, 500)

    h['MuonElectronEvent_MuonPt_SingleMuon'] = ROOT.TH1F( "MuonElectronEvent_MuonPt_SingleMuon" , "MuonElectronEvent_MuonPt_SingleMuon; P_{T} (GeV), Events", 500, 0, 500)
    h['MuonElectronEvent_ElectronPt_SingleMuon'] = ROOT.TH1F( "MuonElectronEvent_ElectronPt_SingleMuon" , "MuonElectronEvent_ElectronPt_SingleMuon; P_{T} (GeV), Events", 500, 0, 500)
    h['MuonElectronEvent_MetPt_SingleMuon'] = ROOT.TH1F( "MuonElectronEvent_MetPt_SingleMuon" , "MuonElectronEvent_MetPt_SingleMuon; P_{T} (GeV), Events", 500, 0, 500)
    h['MuonElectronEvent_dRl_SingleMuon'] = ROOT.TH1F( "MuonElectronEvent_dRl_SingleMuon" , "MuonElectronEvent_dRl_SingleMuon; P_{T} (GeV), Events", 100, 0, 5)
    h['MuonElectronEvent_Mass_SingleMuon'] = ROOT.TH1F( "MuonElectronEvent_Mass_SingleMuon" , "MuonElectronEvent_Mass_SingleMuon; Vis. Mass (GeV), Events", 100, 0, 100)

    h['MuonElectronEvent_MuonPt_MuonEG'] = ROOT.TH1F( "MuonElectronEvent_MuonPt_MuonEG" , "MuonElectronEvent_MuonPt_MuonEG; P_{T} (GeV), Events", 500, 0, 500)
    h['MuonElectronEvent_ElectronPt_MuonEG'] = ROOT.TH1F( "MuonElectronEvent_ElectronPt_MuonEG" , "MuonElectronEvent_ElectronPt_MuonEG; P_{T} (GeV), Events", 500, 0, 500)
    h['MuonElectronEvent_MetPt_MuonEG'] = ROOT.TH1F( "MuonElectronEvent_MetPt_MuonEG" , "MuonElectronEvent_MetPt_MuonEG; P_{T} (GeV), Events", 500, 0, 500)
    h['MuonElectronEvent_dRl_MuonEG'] = ROOT.TH1F( "MuonElectronEvent_dRl_MuonEG" , "MuonElectronEvent_dRl_MuonEG; P_{T} (GeV), Events", 100, 0, 5)
    h['MuonElectronEvent_Mass_MuonEG'] = ROOT.TH1F( "MuonElectronEvent_Mass_MuonEG" , "MuonElectronEvent_Mass_MuonEG; Vis. Mass (GeV), Events", 100, 0, 100)

    h['MuonElectronEvent_MuonPt_Both'] = ROOT.TH1F( "MuonElectronEvent_MuonPt_Both" , "MuonElectronEvent_MuonPt_Both; P_{T} (GeV), Events", 500, 0, 500)
    h['MuonElectronEvent_ElectronPt_Both'] = ROOT.TH1F( "MuonElectronEvent_ElectronPt_Both" , "MuonElectronEvent_ElectronPt_Both; P_{T} (GeV), Events", 500, 0, 500)
    h['MuonElectronEvent_MetPt_Both'] = ROOT.TH1F( "MuonElectronEvent_MetPt_Both" , "MuonElectronEvent_MetPt_Both; P_{T} (GeV), Events", 500, 0, 500)
    h['MuonElectronEvent_dRl_Both'] = ROOT.TH1F( "MuonElectronEvent_dRl_Both" , "MuonElectronEvent_dRl_Both; P_{T} (GeV), Events", 100, 0, 5)
    h['MuonElectronEvent_Mass_Both'] = ROOT.TH1F( "MuonElectronEvent_Mass_Both" , "MuonElectronEvent_Mass_Both; Vis. Mass (GeV), Events", 100, 0, 100)

    h['MuonElectronJetEvent_MuonPt_SingleMuon'] = ROOT.TH1F( "MuonElectronJetEvent_MuonPt_SingleMuon" , "MuonElectronJetEvent_MuonPt_SingleMuon; P_{T} (GeV), Events", 500, 0, 500)
    h['MuonElectronJetEvent_ElectronPt_SingleMuon'] = ROOT.TH1F( "MuonElectronJetEvent_ElectronPt_SingleMuon" , "MuonElectronJetEvent_ElectronPt_SingleMuon; P_{T} (GeV), Events", 500, 0, 500)
    h['MuonElectronJetEvent_JetPt_SingleMuon'] = ROOT.TH1F( "MuonElectronJetEvent_JetPt_SingleMuon" , "MuonElectronJetEvent_JetPt_SingleMuon; P_{T} (GeV), Events", 2000, 0, 2000)
    h['MuonElectronJetEvent_MetPt_SingleMuon'] = ROOT.TH1F( "MuonElectronJetEvent_MetPt_SingleMuon" , "MuonElectronJetEvent_MetPt_SingleMuon; P_{T} (GeV), Events", 500, 0, 500)
    h['MuonElectronJetEvent_dRl_SingleMuon'] = ROOT.TH1F( "MuonElectronJetEvent_dRl_SingleMuon" , "MuonElectronJetEvent_dRl_SingleMuon; P_{T} (GeV), Events", 100, 0, 5)
    h['MuonElectronJetEvent_dRj_SingleMuon'] = ROOT.TH1F( "MuonElectronJetEvent_dRj_SingleMuon" , "MuonElectronJetEvent_dRj_SingleMuon; P_{T} (GeV), Events", 100, 0, 5)
    h['MuonElectronJetEvent_deepjet_SingleMuon'] = ROOT.TH1F( 'MuonElectronJetEvent_deepjet_SingleMuon', 'MuonElectronJetEvent_deepjet_SingleMuon; deepjet; Events', 100, 0, 1)
    h['MuonElectronJetEvent_Mass_SingleMuon'] = ROOT.TH1F( "MuonElectronJetEvent_Mass_SingleMuon" , "MuonElectronJetEvent_Mass_SingleMuon; Vis. Mass (GeV), Events", 100, 0, 100)

    h['MuonElectronJetEvent_MuonPt_MuonEG'] = ROOT.TH1F( "MuonElectronJetEvent_MuonPt_MuonEG" , "MuonElectronJetEvent_MuonPt_MuonEG; P_{T} (GeV), Events", 500, 0, 500)
    h['MuonElectronJetEvent_ElectronPt_MuonEG'] = ROOT.TH1F( "MuonElectronJetEvent_ElectronPt_MuonEG" , "MuonElectronJetEvent_ElectronPt_MuonEG; P_{T} (GeV), Events", 500, 0, 500)
    h['MuonElectronJetEvent_JetPt_MuonEG'] = ROOT.TH1F( "MuonElectronJetEvent_JetPt_MuonEG" , "MuonElectronJetEvent_JetPt_MuonEG; P_{T} (GeV), Events", 2000, 0, 2000)
    h['MuonElectronJetEvent_MetPt_MuonEG'] = ROOT.TH1F( "MuonElectronJetEvent_MetPt_MuonEG" , "MuonElectronJetEvent_MetPt_MuonEG; P_{T} (GeV), Events", 500, 0, 500)
    h['MuonElectronJetEvent_dRl_MuonEG'] = ROOT.TH1F( "MuonElectronJetEvent_dRl_MuonEG" , "MuonElectronJetEvent_dRl_MuonEG; P_{T} (GeV), Events", 100, 0, 5)
    h['MuonElectronJetEvent_dRj_MuonEG'] = ROOT.TH1F( "MuonElectronJetEvent_dRj_MuonEG" , "MuonElectronJetEvent_dRj_MuonEG; P_{T} (GeV), Events", 100, 0, 5)
    h['MuonElectronJetEvent_deepjet_MuonEG'] = ROOT.TH1F( 'MuonElectronJetEvent_deepjet_MuonEG', 'MuonElectronJetEvent_deepjet_MuonEG; deepjet; Events', 100, 0, 1)
    h['MuonElectronJetEvent_Mass_MuonEG'] = ROOT.TH1F( "MuonElectronJetEvent_Mass_MuonEG" , "MuonElectronJetEvent_Mass_MuonEG; Vis. Mass (GeV), Events", 100, 0, 100)

    h['MuonElectronJetEvent_MuonPt_Both'] = ROOT.TH1F( "MuonElectronJetEvent_MuonPt_Both" , "MuonElectronJetEvent_MuonPt_Both; P_{T} (GeV), Events", 500, 0, 500)
    h['MuonElectronJetEvent_ElectronPt_Both'] = ROOT.TH1F( "MuonElectronJetEvent_ElectronPt_Both" , "MuonElectronJetEvent_ElectronPt_Both; P_{T} (GeV), Events", 500, 0, 500)
    h['MuonElectronJetEvent_JetPt_Both'] = ROOT.TH1F( "MuonElectronJetEvent_JetPt_Both" , "MuonElectronJetEvent_JetPt_Both; P_{T} (GeV), Events", 2000, 0, 2000)
    h['MuonElectronJetEvent_MetPt_Both'] = ROOT.TH1F( "MuonElectronJetEvent_MetPt_Both" , "MuonElectronJetEvent_MetPt_Both; P_{T} (GeV), Events", 500, 0, 500)
    h['MuonElectronJetEvent_dRl_Both'] = ROOT.TH1F( "MuonElectronJetEvent_dRl_Both" , "MuonElectronJetEvent_dRl_Both; P_{T} (GeV), Events", 100, 0, 5)
    h['MuonElectronJetEvent_dRj_Both'] = ROOT.TH1F( "MuonElectronJetEvent_dRj_Both" , "MuonElectronJetEvent_dRj_Both; P_{T} (GeV), Events", 100, 0, 5)
    h['MuonElectronJetEvent_deepjet_Both'] = ROOT.TH1F( 'MuonElectronJetEvent_deepjet_Both', 'MuonElectronJetEvent_deepjet_Both; deepjet; Events', 100, 0, 1)
    h['MuonElectronJetEvent_Mass_Both'] = ROOT.TH1F( "MuonElectronJetEvent_Mass_Both" , "MuonElectronJetEvent_Mass_Both; Vis. Mass (GeV), Events", 100, 0, 100)

    h['MuonElectronJethighMetEvent_MuonPt_SingleMuon'] = ROOT.TH1F( "MuonElectronJethighMetEvent_MuonPt_SingleMuon" , "MuonElectronJethighMetEvent_MuonPt_SingleMuon; P_{T} (GeV), Events", 500, 0, 500)
    h['MuonElectronJethighMetEvent_ElectronPt_SingleMuon'] = ROOT.TH1F( "MuonElectronJethighMetEvent_ElectronPt_SingleMuon" , "MuonElectronJethighMetEvent_ElectronPt_SingleMuon;P_{T} (GeV), Events", 500, 0, 500)
    h['MuonElectronJethighMetEvent_JetPt_SingleMuon'] = ROOT.TH1F( "MuonElectronJethighMetEvent_JetPt_SingleMuon" , "MuonElectronJethighMetEvent_JetPt_SingleMuon; P_{T} (GeV), Events", 2000, 0, 2000)
    h['MuonElectronJethighMetEvent_dRl_SingleMuon'] = ROOT.TH1F( "MuonElectronJethighMetEvent_dRl_SingleMuon" , "MuonElectronJethighMetEvent_dRl_SingleMuon; P_{T} (GeV), Events", 100, 0, 5)
    h['MuonElectronJethighMetEvent_dRj_SingleMuon'] = ROOT.TH1F( "MuonElectronJethighMetEvent_dRj_SingleMuon" , "MuonElectronJethighMetEvent_dRj_SingleMuon; P_{T} (GeV), Events", 100, 0, 5)
    h['MuonElectronJethighMetEvent_deepjet_SingleMuon'] = ROOT.TH1F( 'MuonElectronJethighMetEvent_deepjet_SingleMuon', 'MuonElectronJethighMetEvent_deepjet_SingleMuon; deepjet; Events', 100, 0, 1)
    h['MuonElectronJethighMetEvent_Mass_SingleMuon'] = ROOT.TH1F( "MuonElectronJethighMetEvent_Mass_SingleMuon" , "MuonElectronJethighMetEvent_Mass_SingleMuon; Vis. Mass (GeV), Events", 100, 0, 100)

    h['MuonElectronJethighMetEvent_MuonPt_MuonEG'] = ROOT.TH1F( "MuonElectronJethighMetEvent_MuonPt_MuonEG" , "MuonElectronJethighMetEvent_MuonPt_MuonEG; P_{T} (GeV), Events", 500, 0, 500)
    h['MuonElectronJethighMetEvent_ElectronPt_MuonEG'] = ROOT.TH1F( "MuonElectronJethighMetEvent_ElectronPt_MuonEG" , "MuonElectronJethighMetEvent_ElectronPt_MuonEG; P_{T} (GeV), Events", 500, 0, 500)
    h['MuonElectronJethighMetEvent_JetPt_MuonEG'] = ROOT.TH1F( "MuonElectronJethighMetEvent_JetPt_MuonEG" , "MuonElectronJethighMetEvent_JetPt_MuonEG; P_{T} (GeV), Events", 2000, 0, 2000)
    h['MuonElectronJethighMetEvent_dRl_MuonEG'] = ROOT.TH1F( "MuonElectronJethighMetEvent_dRl_MuonEG" , "MuonElectronJethighMetEvent_dRl_MuonEG; P_{T} (GeV), Events", 100, 0, 5)
    h['MuonElectronJethighMetEvent_dRj_MuonEG'] = ROOT.TH1F( "MuonElectronJethighMetEvent_dRj_MuonEG" , "MuonElectronJethighMetEvent_dRj_MuonEG; P_{T} (GeV), Events", 100, 0, 5)
    h['MuonElectronJethighMetEvent_deepjet_MuonEG'] = ROOT.TH1F( 'MuonElectronJethighMetEvent_deepjet_MuonEG', 'MuonElectronJethighMetEvent_deepjet_MuonEG; deepjet; Events', 100, 0, 1)
    h['MuonElectronJethighMetEvent_Mass_MuonEG'] = ROOT.TH1F( "MuonElectronJethighMetEvent_Mass_MuonEG" , "MuonElectronJethighMetEvent_Mass_MuonEG; Vis. Mass (GeV), Events", 100, 0, 100)

    h['MuonElectronJethighMetEvent_MuonPt_Both'] = ROOT.TH1F( "MuonElectronJethighMetEvent_MuonPt_Both" , "MuonElectronJethighMetEvent_MuonPt_Both; P_{T} (GeV), Events", 500, 0,500)
    h['MuonElectronJethighMetEvent_ElectronPt_Both'] = ROOT.TH1F( "MuonElectronJethighMetEvent_ElectronPt_Both" , "MuonElectronJethighMetEvent_ElectronPt_Both; P_{T} (GeV), Events", 500, 0, 500)
    h['MuonElectronJethighMetEvent_JetPt_Both'] = ROOT.TH1F( "MuonElectronJethighMetEvent_JetPt_Both" , "MuonElectronJethighMetEvent_JetPt_Both; P_{T} (GeV), Events", 2000, 0, 2000)
    h['MuonElectronJethighMetEvent_dRl_Both'] = ROOT.TH1F( "MuonElectronJethighMetEvent_dRl_Both" , "MuonElectronJethighMetEvent_dRl_Both; P_{T} (GeV), Events", 100, 0, 5)
    h['MuonElectronJethighMetEvent_dRj_Both'] = ROOT.TH1F( "MuonElectronJethighMetEvent_dRj_Both" , "MuonElectronJethighMetEvent_dRj_Both; P_{T} (GeV), Events", 100, 0, 5)
    h['MuonElectronJethighMetEvent_deepjet_Both'] = ROOT.TH1F( 'MuonElectronJethighMetEvent_deepjet_Both', 'MuonElectronJethighMetEvent_deepjet_Both; deepjet; Events', 100, 0, 1)
    h['MuonElectronJethighMetEvent_Mass_Both'] = ROOT.TH1F( "MuonElectronJethighMetEvent_Mass_Both" , "MuonElectronJethighMetEvent_Mass_Both; Vis. Mass (GeV), Events", 100, 0, 100)


hist_regions = []

trigger_regions = [
    "_lowMET_SingleMuon",
    "_bjet_SingleMuon",
    "_lowMET_MuonEG",
    "_bjet_MuonEG",
    "_lowMET_Both_SingleMuon",
    "_bjet_Both_SingleMuon",
    "_lowMET_Both_MuonEG",
    "_bjet_Both_MuonEG"
]

Regions = [
    "_baseline",
    "_baseline_bjet0",
    "_baseline_bjet1",
    "_baseline_bjet2",
    "_baseline_bjet3",
    "_dRcut",
    "_SR",
    "_dRcut_dPhicut_lowMET",
    "_dRcut_bjet",
    "_dRcut_highMET",
    "_dRcut_highMET_bjet0",
    "_dRcut_highMET_bjet",
    "_dRcut_highMET_bjet1",
    "_dRcut_highMET_bjet2",
    "_dRcut_highMET_bjet3",
    "_dRcut_bjet0",
    "_dRcut_bjet1",
    "_dRcut_bjet2",
    "_dRcut_bjet3",
    "_dRcut_lowMET",
    "_dRcut_lowMET_bjet0",
    "_dRcut_lowMET_bjet",
    "_dRcut_lowMET_bjet1",
    "_dRcut_lowMET_bjet2",
    "_dRcut_lowMET_bjet3",
    "_dRcut_highMET_dPhicut_bjet",
    "_dRcut_highMET_dPhicut_bjet1",
    "_dRcut_highMET_dPhicut_bjet2",
    "_dRcut_highMET_dPhicut_bjet3",
    "_dRcut_lowMET_dPhicut_bjet",
    "_lowMET_bjet0",
    "_lowMET_bjet1",
    "_lowMET_bjet2",
    "_lowMET_bjet3",
    "_highMET_bjet0",
    "_highMET_bjet1",
    "_highMET_bjet2",
    "_highMET_bjet3",
    "_lowMET",
    "_highMET",
    "_bjet",
    "_SR_dPhicut",
    "_SR_mTcut"
]

Bjet_Regions = [
    "_2bjet",
    "_lowMET_2bjet",
    "_highMET_2bjet"
]

bjet_hist_regions = []

if triggerStudy == True:
    Regions.extend(trigger_regions)

b_eff_regions = [
    "_baseline",
    "_dRcut",
    "_dRcut_lowMET",
    "_dRcut_highMET",
    "_lowMET",
    "_highMET",
]


isoRegions = ["hEMu",
              "hEMu_nonIsoMu",
              "hEMu_nonIsoE",
              "hEMu_nonIso",
              "hEMu_checkMuMu",
]

isoPlotRegions = []

for rb in Bjet_Regions:
    bjet_hist_regions.append("hMuMu"+rb)

for isoR in isoRegions :
    for rb in Bjet_Regions:
        bjet_hist_regions.append(isoR+rb)
    for re in b_eff_regions:
        define_b_eff_histogram(isoR+re)
    for r in Regions:
        hist_regions.append(isoR+r)
        if r == "_SR" or r == "_SR_dPhicut" or r == "_dRcut_highMET_bjet" or r == "_dRcut_highMET_dPhicut_bjet" :
            isoPlotRegions.append(isoR+r)

kin_regions = [
    "_baseline",
    "_ePtcut",
    "_lPtcut",
    "_jetPtcut"
]


mumu_regions = [
    "hMuMu_Baseline",
    "hMuMu_dRcut",
    "hMuMu_dRcut_highMET",
    "hMuMu_SR_dRcut_highMET_dPhicut"
]

hist_regions.extend(mumu_regions)

for rb in bjet_hist_regions:
    define_bjet_histogram(rb)

for r in hist_regions:
    define_event_histogram(r)

for rk in kin_regions:
    define_event_histogram("kin"+rk)

define_general_histogram()
cut_flow_hists()

for key in h.keys():
    h[key].Sumw2()



def get_TLorentzVector(l1, l2, js, met_pt, met_phi):

    vl1 = ROOT.TLorentzVector()
    vl1.SetPtEtaPhiM(l1.pt, l1.eta, l1.phi, l1.mass)

    vl2 = ROOT.TLorentzVector()
    vl2.SetPtEtaPhiM(l2.pt, l2.eta, l2.phi, l2.mass)

    vj = ROOT.TLorentzVector()
    vj.SetPtEtaPhiM(js.pt, js.eta, js.phi, js.mass)

    vm = ROOT.TLorentzVector()
    vm.SetPtEtaPhiM(met_pt, 0, met_phi, 0)

    return vl1, vl2, vj, vm


def pass_deltaR(l1,l2,j, channel):

    if channel == "muTau" or channel == "eTau":
        if l1.DeltaR(l2) < event_cut["dRl"] and j.DeltaR(l1) > event_cut["dRlj"] and j.DeltaR(l2) > event_cut["dRlj"] and l1.DeltaR(l2) > event_cut["dRltau"]:
            return 1
        else:
            return -9999
    else:
        if l1.DeltaR(l2) < event_cut["dRl"] and j.DeltaR(l1) > event_cut["dRlj"] and j.DeltaR(l2) > event_cut["dRlj"]:
            return 1
        else:
            return -9999


def pass_baseline(l1, l2, j):

    if j.Pt() > event_cut['jetPt'] and (l1+l2).M() > event_cut["mass"]:
        return 1
    else:
        return -9999


def MuMu_Channel(mu, js, met_pt, met_phi, plot=False):

    isMuMu = 0
    h['hMuMu_Events'].Fill(1, genweight)

    mu1, mu2, j, m = get_TLorentzVector(mu[0], mu[1], js[0], met_pt, met_phi)

    if pass_baseline(mu1, mu2, j):

        isJetHTEvent = 0
        isSingleMuonEvent = 0

        if ( j.Pt() > 500 and isHT == 1 ) : isJetHTEvent = 1

        if ( mu1.Pt() > 50 and isMu == 1 ) or ( mu1.Pt() > 27 and isIsoMu == 1 ) : isSingleMuonEvent = 1

        if ( isData == 0 and isSingleMuonEvent == 1 ) or \
           ( isData == 1 and isSingleMuonDataset == 1 and isSingleMuonEvent == 1 ) :

            if plot == True:

#                check_Nbjet(mu1, mu2, j, m, js, "hMuMu")            
                plot_event_hist("hMuMu_Baseline", mu1, mu2, j, m)

            if pass_deltaR(mu1, mu2, j, "MuMu"):

                if plot==True: plot_event_hist('hMuMu_dRcut', mu1, mu2, j, m)

                if m.Pt() > event_cut['metcut']:
                    if plot==True: plot_event_hist('hMuMu_dRcut_highMET', mu1, mu2, j, m)

                    if abs(m.DeltaPhi(mu1)) < 1 and abs(m.DeltaPhi(j)) > 2:
                        if plot==True: plot_event_hist('hMuMu_SR_dRcut_highMET_dPhicut', mu1, mu2, j, m)

                        isMuMu = 1

    return isMuMu


def EMu_Channel(ele, mu_emu, js, met_pt, met_phi, isolation = "hEMu"):

    isEMu = 0
    h[isolation+'_Events'].Fill(1, genweight)

    e, mu, j, m = get_TLorentzVector(ele[0], mu_emu[0], js[0], met_pt, met_phi)

    trigger = [0,0]

    if ( mu.Pt() > 50 and isMu == 1 ) or ( mu.Pt() > 27 and isIsoMu == 1 ) : trigger[0] = 1
    
    if ( ( ( mu.Pt() > 8 and e.Pt() > 23 ) or ( mu.Pt() > 23 and e.Pt() > 12 ) ) and isMuonEG == 1 ) : trigger[1] = 1

    if triggerStudy == True :

        trigger_study(e, mu, j , m, trigger)

    if ( isData == 0 and ( trigger[0] == 1 or trigger[1] == 1 ) ) or \
       ( isData == 1 and ( isSingleMuonDataset == 1 and trigger[0] == 1 ) ) or \
       ( isData == 1 and ( isMuonEGDataset == 1 and ( trigger[0] == 0 and trigger[1] == 1 ) ) ) :

        if pass_baseline(e, mu, j) == 1 :
            plot_event_hist(isolation+"_baseline", e, mu, j, m)

            check_Nbjet(e, mu, j, m, js, isolation)

            if len(s_b) == 0 : plot_event_hist(isolation+"_baseline_bjet0", e, mu, j, m)
            if len(s_b) == 1 : plot_event_hist(isolation+"_baseline_bjet1", e, mu, j, m)
            if len(s_b) == 2 : plot_event_hist(isolation+"_baseline_bjet2", e, mu, j, m)
            if len(s_b) > 1 : plot_event_hist(isolation+"_baseline_bjet3", e, mu, j, m)

            if len(s_b) == 0 : h[isolation+'_Events'].Fill(2, genweight)

            if m.Pt() < event_cut['metcut'] : #lowMET
                plot_event_hist(isolation+"_lowMET", e, mu, j, m)

                if len(s_b) == 0 : plot_event_hist(isolation+"_lowMET_bjet0", e, mu, j, m)
                if len(s_b) == 1 : plot_event_hist(isolation+"_lowMET_bjet1", e, mu, j, m)
                if len(s_b) == 2 : plot_event_hist(isolation+"_lowMET_bjet2", e, mu, j, m)
                if len(s_b) > 1 : plot_event_hist(isolation+"_lowMET_bjet3", e, mu, j, m)

            if m.Pt() > event_cut['metcut'] : #highMET                                                                                                                     
                plot_event_hist(isolation+"_highMET", e, mu, j, m)

                if len(s_b) == 0 : plot_event_hist(isolation+"_highMET_bjet0", e, mu, j, m)
                if len(s_b) == 1 : plot_event_hist(isolation+"_highMET_bjet1", e, mu, j, m)
                if len(s_b) == 2 : plot_event_hist(isolation+"_highMET_bjet2", e, mu, j, m)
                if len(s_b) > 1 : plot_event_hist(isolation+"_highMET_bjet3", e, mu, j, m)

            if len(s_b) > 0 : #non-zero b-jet
                plot_event_hist(isolation+"_bjet", e, mu, j, m)

            if pass_deltaR(e, mu, j, "EMu") == 1: #dRcut

                plot_event_hist(isolation+"_dRcut", e, mu, j ,m)
                find_b_efficiency(isolation+"_dRcut", js)

                if len(s_b) == 0 :
                    plot_event_hist(isolation+"_dRcut_bjet0", e, mu, j ,m)
                
                if m.Pt() > event_cut['metcut'] and ( isData == 0 or isNonIso == True ) : 
                    h[isolation+'_Events'].Fill(3, genweight)

                    plot_event_hist(isolation+"_dRcut_highMET", e, mu, j ,m) #dRcut highMET

                    if len(s_b) == 0 : #SR
                        h[isolation+'_Events'].Fill(4, genweight)
                        plot_event_hist(isolation+"_SR", e, mu, j, m)

                        h[isolation+"_SR_muIso"].Fill(mu_emu[0].iso, genweight)
                            
                        if np.cos(m.DeltaPhi(mu)) > 0.8 : #SR dPhicut
                            h[isolation+'_Events'].Fill(5, genweight)
                            plot_event_hist(isolation+"_SR_dPhicut", e, mu, j, m)

                            h[isolation+"_SR_dPhicut_muIso"].Fill(mu_emu[0].iso, genweight)

                        if Mt((e+mu),m) < 40 :
                            h[isolation+'_Events'].Fill(6, genweight)
                            plot_event_hist(isolation+"_SR_mTcut", e, mu, j, m)
                    
                if m.Pt() < event_cut['metcut'] : #dRcut lowMET
                    plot_event_hist(isolation+"_dRcut_lowMET", e, mu, j, m)

                    if len(s_b) == 0 :
                        plot_event_hist(isolation+"_dRcut_lowMET_bjet0", e, mu, j, m)

                    if np.cos(m.DeltaPhi(mu)) > 0.8 :
                        plot_event_hist(isolation+"_dRcut_dPhicut_lowMET", e, mu, j, m)

                if len(s_b) > 0 : #dRcut non-zero b-jet
                    plot_event_hist(isolation+"_dRcut_bjet", e, mu, j, m)

                    if len(s_b) == 1 : plot_event_hist(isolation+"_dRcut_bjet1", e, mu, j, m)
                    if len(s_b) == 2 : plot_event_hist(isolation+"_dRcut_bjet2", e, mu, j, m)
                    if len(s_b) > 1 : plot_event_hist(isolation+"_dRcut_bjet3", e, mu, j, m)

                    if m.Pt() > event_cut['metcut'] : #high MET non-zero b-jet
                        plot_event_hist(isolation+"_dRcut_highMET_bjet", e, mu, j, m)

                        if len(s_b) == 1 : plot_event_hist(isolation+"_dRcut_highMET_bjet1", e, mu, j, m)
                        if len(s_b) == 2 : plot_event_hist(isolation+"_dRcut_highMET_bjet2", e, mu, j, m)
                        if len(s_b) > 1 : plot_event_hist(isolation+"_dRcut_highMET_bjet3", e, mu, j, m)

                        h[isolation+"_dRcut_highMET_bjet_muIso"].Fill(mu_emu[0].iso, genweight)

                        if np.cos(m.DeltaPhi(mu)) > 0.8 :
                            plot_event_hist(isolation+"_dRcut_highMET_dPhicut_bjet", e, mu, j, m)

                            if len(s_b) == 1 : plot_event_hist(isolation+"_dRcut_highMET_dPhicut_bjet1", e, mu, j, m)
                            if len(s_b) == 2 : plot_event_hist(isolation+"_dRcut_highMET_dPhicut_bjet2", e, mu, j, m)
                            if len(s_b) > 1 : plot_event_hist(isolation+"_dRcut_highMET_dPhicut_bjet3", e, mu, j, m)

                            h[isolation+"_dRcut_highMET_dPhicut_bjet_muIso"].Fill(mu_emu[0].iso, genweight)

                    if m.Pt() < event_cut['metcut'] : #lowMET non-zero b-jet
                        plot_event_hist(isolation+"_dRcut_lowMET_bjet", e, mu, j, m)

                        if len(s_b) == 1 : plot_event_hist(isolation+"_dRcut_lowMET_bjet1", e, mu, j, m)
                        if len(s_b) == 2 : plot_event_hist(isolation+"_dRcut_lowMET_bjet2", e, mu, j, m)
                        if len(s_b) > 1 : plot_event_hist(isolation+"_dRcut_lowMET_bjet3", e, mu, j, m)

                        if np.cos(m.DeltaPhi(mu)) > 0.8 :
                            plot_event_hist(isolation+"_dRcut_lowMET_dPhicut_bjet", e, mu, j, m)



def find_b_efficiency(isolation, jets):

    if bEffStudy == True :

        for i in range(len(jets)):
            if jets[i].jetflavour == 5 :
                h[isolation+'_BFlavour_JetPt'].Fill(jets[i].pt)
                h[isolation+'_BFlavour_JetEta'].Fill(abs(jets[i].eta))
                h[isolation+'_BFlavour_JetPtEta'].Fill(jets[i].pt, abs(jets[i].eta))
                if jets[i].deepjet > 0.7476 :
                    h[isolation+'_BFlavour_BTagged_JetPt'].Fill(jets[i].pt)
                    h[isolation+'_BFlavour_BTagged_JetEta'].Fill(abs(jets[i].eta))
                    h[isolation+'_BFlavour_BTagged_JetPtEta'].Fill(jets[i].pt, abs(jets[i].eta))

            if jets[i].jetflavour == 4 :
                h[isolation+'_CFlavour_JetPt'].Fill(jets[i].pt)
                h[isolation+'_CFlavour_JetEta'].Fill(abs(jets[i].eta))
                h[isolation+'_CFlavour_JetPtEta'].Fill(jets[i].pt, abs(jets[i].eta))
                if jets[i].deepjet > 0.7476:
                    h[isolation+'_CFlavour_BTagged_JetPt'].Fill(jets[i].pt)
                    h[isolation+'_CFlavour_BTagged_JetEta'].Fill(abs(jets[i].eta))
                    h[isolation+'_CFlavour_BTagged_JetPtEta'].Fill(jets[i].pt, abs(jets[i].eta))

            if jets[i].jetflavour == 0 :
                h[isolation+'_LFlavour_JetPt'].Fill(jets[i].pt)
                h[isolation+'_LFlavour_JetEta'].Fill(abs(jets[i].eta))
                h[isolation+'_LFlavour_JetPtEta'].Fill(jets[i].pt, abs(jets[i].eta))
                if jets[i].deepjet > 0.7476:
                    h[isolation+'_LFlavour_BTagged_JetPt'].Fill(jets[i].pt)
                    h[isolation+'_LFlavour_BTagged_JetEta'].Fill(abs(jets[i].eta))
                    h[isolation+'_LFlavour_BTagged_JetPtEta'].Fill(jets[i].pt, abs(jets[i].eta))


def plot_for_bjet(region,bj1, bj2, e, mu):

    h[region+'_dRbjet'].Fill(bj1.DeltaR(bj2), genweight)
    h[region+'_dRebjet1'].Fill(bj1.DeltaR(e), genweight)
    h[region+'_dRebjet2'].Fill(bj2.DeltaR(e), genweight)
    h[region+'_dRmubjet1'].Fill(bj1.DeltaR(mu), genweight)
    h[region+'_dRmubjet2'].Fill(bj2.DeltaR(mu), genweight)

def plot_pt_eta_deepjet(name, js):

    for i in range(len(js)):
        if js[i].deepjet > 0.0532 and js[i].deepjet < 0.3040:
            h[name+'_Loosebjet_Pt'].Fill(js[i].pt, genweight)
            h[name+'_Loosebjet_Eta'].Fill(js[i].eta, genweight)
        if js[i].deepjet > 0.3040 and js[i].deepjet < 0.7476:
            h[name+'_Mediumbjet_Pt'].Fill(js[i].pt, genweight)
            h[name+'_Mediumbjet_Eta'].Fill(js[i].eta, genweight)
        if js[i].deepjet > 0.7476 :
            h[name+'_Tightbjet_Pt'].Fill(js[i].pt, genweight)
            h[name+'_Tightbjet_Eta'].Fill(js[i].eta, genweight)


def check_Nbjet(e, mu, jv, m, js, isolation = "hEMu"):

    h[isolation+'_Nbjet'].Fill(len(s_b), genweight)
    find_b_efficiency(isolation+"_baseline", js)
    plot_pt_eta_deepjet(isolation+"_baseline", js)

    for i in range(len(js)):
        h[isolation+'_deepjet'].Fill(js[i].deepjet, genweight)

    if len(s_b) > 0 :
        bj = ROOT.TLorentzVector()
        bj.SetPtEtaPhiM(s_b[0].pt, s_b[0].eta, s_b[0].phi, s_b[0].mass)

        h[isolation+'_deepjet_dRlbj'].Fill(s_b[0].deepjet, mu.DeltaR(bj), genweight)
        
    if len(s_b) == 2 :

        bj1 = ROOT.TLorentzVector()
        bj1.SetPtEtaPhiM(s_b[0].pt, s_b[0].eta, s_b[0].phi, s_b[0].mass)

        bj2 = ROOT.TLorentzVector()
        bj2.SetPtEtaPhiM(s_b[1].pt, s_b[1].eta, s_b[1].phi, s_b[1].mass)

        plot_for_bjet(isolation+'_2bjet', bj1, bj2, e, mu)

        for i in range(len(s_b)):
            h[isolation+'_2bjet_deepjet'].Fill(s_b[i].deepjet, genweight)

    if e.DeltaR(mu) < 0.8 and jv.DeltaR(e+mu) > 0.8 :

        for i in range(len(js)):
            h[isolation+'_loosedRlcut_deepjet'].Fill(js[i].deepjet, genweight)

    if jv.DeltaR(e+mu) > 0.8 :

        for i in range(len(js)):
            h[isolation+'_nodRlcut_deepjet'].Fill(js[i].deepjet, genweight)
    
    if pass_deltaR(e, mu, jv, "EMu") == 1:

        h[isolation+'_dRcut_Nbjet'].Fill(len(s_b), genweight)
        find_b_efficiency(isolation+"_dRcut", js)
        plot_pt_eta_deepjet(isolation+"_dRcut", js)

        for i in range(len(js)):
            h[isolation+'_dRcut_deepjet'].Fill(js[i].deepjet, genweight)

        if len(s_b) > 0 :
            h[isolation+'_dRcut_deepjet_dRlbj'].Fill(s_b[0].deepjet, mu.DeltaR(bj), genweight)
        
        if m.Pt() < event_cut['metcut'] :

            h[isolation+'_dRcut_lowMET_Nbjet'].Fill(len(s_b), genweight)
            find_b_efficiency(isolation+"_dRcut_lowMET", js)
            plot_pt_eta_deepjet(isolation+"_dRcut_lowMET", js)

            for i in range(len(js)):
                h[isolation+'_dRcut_lowMET_deepjet'].Fill(js[i].deepjet, genweight)

            if len(s_b) > 0:
                h[isolation+'_dRcut_lowMET_deepjet_dRlbj'].Fill(s_b[0].deepjet, mu.DeltaR(bj), genweight)

            if np.cos(m.DeltaPhi(mu)) > 0.8 :

                h[isolation+'_dRcut_lowMET_dPhicut_Nbjet'].Fill(len(s_b), genweight)
                for i in range(len(js)):
                    h[isolation+'_dRcut_lowMET_dPhicut_deepjet'].Fill(js[i].deepjet, genweight)

        if m.Pt() > event_cut['metcut'] and isData == 0:

            h[isolation+'_dRcut_highMET_Nbjet'].Fill(len(s_b), genweight)
            find_b_efficiency(isolation+"_dRcut_highMET", js)
            plot_pt_eta_deepjet(isolation+"_dRcut_highMET", js)

            for i in range(len(js)):
                h[isolation+'_dRcut_highMET_deepjet'].Fill(js[i].deepjet, genweight)

            if len(s_b) > 0:
                h[isolation+'_dRcut_highMET_deepjet_dRlbj'].Fill(s_b[0].deepjet, mu.DeltaR(bj), genweight)

    if m.Pt() < event_cut['metcut'] :

        if len(s_b) == 2 :

            bj1 = ROOT.TLorentzVector()
            bj1.SetPtEtaPhiM(s_b[0].pt, s_b[0].eta, s_b[0].phi, s_b[0].mass)

            bj2 = ROOT.TLorentzVector()
            bj2.SetPtEtaPhiM(s_b[1].pt, s_b[1].eta, s_b[1].phi, s_b[1].mass)

            plot_for_bjet(isolation+'_lowMET_2bjet', bj1, bj2, e, mu)

            for i in range(len(s_b)):
                h[isolation+'_lowMET_2bjet_deepjet'].Fill(s_b[i].deepjet, genweight)

        h[isolation+'_lowMET_Nbjet'].Fill(len(s_b), genweight)
        find_b_efficiency(isolation+"_lowMET", js)
        plot_pt_eta_deepjet(isolation+"_lowMET", js)

        for i in range(len(js)):
            h[isolation+'_lowMET_deepjet'].Fill(js[i].deepjet, genweight)

        if len(s_b) > 0:
            h[isolation+'_lowMET_deepjet_dRlbj'].Fill(s_b[0].deepjet, mu.DeltaR(bj), genweight)

    if m.Pt() > event_cut['metcut'] :

        if len(s_b) == 2 :

            bj1 = ROOT.TLorentzVector()
            bj1.SetPtEtaPhiM(s_b[0].pt, s_b[0].eta, s_b[0].phi, s_b[0].mass)

            bj2 = ROOT.TLorentzVector()
            bj2.SetPtEtaPhiM(s_b[1].pt, s_b[1].eta, s_b[1].phi, s_b[1].mass)

            plot_for_bjet(isolation+'_highMET_2bjet', bj1, bj2, e, mu)

            for i in range(len(s_b)):
                h[isolation+'_highMET_2bjet_deepjet'].Fill(s_b[i].deepjet, genweight)

        h[isolation+'_highMET_Nbjet'].Fill(len(s_b), genweight)
        find_b_efficiency(isolation+"_highMET", js)
        plot_pt_eta_deepjet(isolation+"_highMET", js)

        for i in range(len(js)):
            h[isolation+'_highMET_deepjet'].Fill(js[i].deepjet, genweight)

        if len(s_b) > 0:
            h[isolation+'_highMET_deepjet_dRlbj'].Fill(s_b[0].deepjet, mu.DeltaR(bj), genweight)


def trigger_study(e, mu, j , m, trigger):

    if trigger[0] == 1 :
        if pass_baseline(e, mu, j) == 1 :
            plot_event_hist(isolation+"_Baseline_SingleMuon", e, mu, j ,m)

    if trigger[1] == 1 :
        if pass_baseline(e, mu, j) == 1 :
            plot_event_hist(isolation+"_Baseline_MuonEG", e, mu, j ,m)

    if ( isData == 0 and ( trigger[0] == 1 and trigger[1] == 0 ) ) or \
       ( isData == 1 and ( isSingleMuonDataset == 1 and trigger[0] == 1 and trigger[1] == 0 ) ) :
        if pass_baseline(e, mu, j) == 1 :
            if m.Pt() < event_cut['metcut'] : #lowMET SingleMuon                                                                                          
                plot_event_hist(isolation+"_lowMET_SingleMuon", e, mu, j, m)
            if len(s_b) > 0 : #bjet SingleMuon
                plot_event_hist(isolation+"_bjet_SingleMuon", e, mu, j, m)

    if ( isData == 0 and ( trigger[0] == 0 and trigger[1] == 1 ) ) or \
       ( isData == 1 and ( isMuonEGDataset == 1 and trigger[0] == 0 and trigger[1] == 1 ) ) :
        if pass_baseline(e, mu, j) == 1 :
            if m.Pt() < event_cut['metcut'] :#lowMET MuonEG                                                                                                            
                plot_event_hist(isolation+"_lowMET_MuonEG", e, mu, j, m)
            if len(s_b) > 0 : #bjet MuonEG
                plot_event_hist(isolation+"_bjet_MuonEG", e, mu, j, m)

    if ( isData == 0 and ( trigger[0] == 1 and trigger[1] == 1 ) ) or \
       ( isData == 1 and ( isSingleMuonDataset == 1 and trigger[0] == 1 and trigger[1] == 1 ) ):
        if pass_baseline(e, mu, j) == 1 :
            if m.Pt() < event_cut['metcut'] :#lowMET Both - SingleMuon                                   
                plot_event_hist(isolation+"_lowMET_Both_SingleMuon", e, mu, j, m)
            if len(s_b) > 0 : #bjet Both - SingleMuon
                plot_event_hist(isolation+"_bjet_Both_SingleMuon", e, mu, j, m)

    if ( isData == 0 and ( trigger[0] == 1 and trigger[1] == 1 ) ) or \
       ( isData == 1 and ( isMuonEGDataset == 1 and trigger[0] == 1 and trigger[1] == 1 ) ) :
        if pass_baseline(e, mu, j) == 1 :
            if m.Pt() < event_cut['metcut'] :#lowMET Both - MuonEG
                plot_event_hist(isolation+"_lowMET_Both_MuonEG", e, mu, j, m)
            if len(s_b) > 0 : #bjet Both - MuonEG                                                                                                         
                plot_event_hist(isolation+"_bjet_Both_MuonEG", e, mu, j, m) 



def kin_selection(muons, electrons, jets, isolation="kin"):

    trigger = [0,0]

    e, mu, j1, m = get_TLorentzVector(electrons[0], muons[0], jets[0], met_pt, met_phi)

    j2 = ROOT.TLorentzVector()
    j2.SetPtEtaPhiM(jets[1].pt, jets[1].eta, jets[1].phi, jets[1].mass)

    if ( ( ( mu.Pt() > 8 and e.Pt() > 23 ) or ( mu.Pt() > 23 and e.Pt() > 12 ) ) and isMuonEG == 1 ) : trigger[1] = 1

    if ( mu.Pt() > 50 and isMu == 1 ) or ( mu.Pt() > 27 and isIsoMu == 1 ) : trigger[0] = 1

    if ( isData == 0 and ( trigger[0] == 1 or trigger[1] == 1 ) ) or \
       ( isData == 1 and ( isSingleMuonDataset == 1 and trigger[0] == 1 ) ) or \
       ( isData == 1 and ( isMuonEGDataset == 1 and ( trigger[0] == 0 and trigger[1] == 1 ) ) ) :

        plot_event_hist(isolation+"_baseline", e, mu, j1 ,m)
        for i in range(len(jets)):
            h['kin_baseline_deepjet'].Fill(jets[i].deepjet, genweight)

        if e.Pt() > 25.0 :
            plot_event_hist(isolation+"_ePtcut", e, mu, j1 ,m)
            for i in range(len(jets)):
                h['kin_ePtcut_deepjet'].Fill(jets[i].deepjet, genweight)

            if mu.Pt() > 25.0 :
                plot_event_hist(isolation+"_lPtcut", e, mu, j1 ,m)
                for i in range(len(jets)):
                    h['kin_lPtcut_deepjet'].Fill(jets[i].deepjet, genweight)

                if j1.Pt() > 30.0 and j2.Pt() > 30.0 :
                    plot_event_hist(isolation+"_jetPtcut", e, mu, j1 ,m)
                    for i in range(len(jets)):
                        h['kin_jetPtcut_deepjet'].Fill(jets[i].deepjet, genweight)


def check_cut_flow(muons, electrons, jets):

    trigger = [0,0]

    mu = ROOT.TLorentzVector()
    mu.SetPtEtaPhiM(muons[0].pt, muons[0].eta, muons[0].phi, muons[0].mass)

    if ( mu.Pt() > 50 and isMu == 1 ) or ( mu.Pt() > 27 and isIsoMu == 1 ) : trigger[0] = 1

    if( isData == 0 and ( trigger[0] == 1 ) ) or \
      ( isData == 1 and ( isSingleMuonDataset == 1 and trigger[0] == 1 ) ) :

        h['MuonEvent_MuonPt_SingleMuon'].Fill(mu.Pt(), genweight)
        h['MuonEvent_MetPt_SingleMuon'].Fill(met_pt, genweight)

        if len(electrons) > 0 and muons[0].charge*electrons[0].charge < 0 :

            e = ROOT.TLorentzVector()
            e.SetPtEtaPhiM(electrons[0].pt, electrons[0].eta, electrons[0].phi, electrons[0].mass)

            if (mu+e).M() > 5.0 :

                h['MuonElectronEvent_MuonPt_SingleMuon'].Fill(mu.Pt(), genweight)
                h['MuonElectronEvent_ElectronPt_SingleMuon'].Fill(e.Pt(), genweight)
                h['MuonElectronEvent_MetPt_SingleMuon'].Fill(met_pt, genweight)
                h['MuonElectronEvent_dRl_SingleMuon'].Fill(mu.DeltaR(e), genweight)
                h['MuonElectronEvent_Mass_SingleMuon'].Fill((mu+e).M(), genweight)

            if len(jets) > 0 :

                j = ROOT.TLorentzVector()
                j.SetPtEtaPhiM(jets[0].pt, jets[0].eta, jets[0].phi, jets[0].mass)

                if (mu+e).M() > 5.0 and j.Pt() > event_cut['jetPt'] :

                    h['MuonElectronJetEvent_MuonPt_SingleMuon'].Fill(mu.Pt(), genweight)
                    h['MuonElectronJetEvent_ElectronPt_SingleMuon'].Fill(e.Pt(), genweight)
                    h['MuonElectronJetEvent_JetPt_SingleMuon'].Fill(j.Pt(), genweight)
                    h['MuonElectronJetEvent_MetPt_SingleMuon'].Fill(met_pt, genweight)
                    h['MuonElectronJetEvent_dRl_SingleMuon'].Fill(mu.DeltaR(e), genweight)
                    h['MuonElectronJetEvent_dRj_SingleMuon'].Fill(j.DeltaR(e+mu), genweight)
                    h['MuonElectronJetEvent_Mass_SingleMuon'].Fill((mu+e).M(), genweight)
                    for i in range(len(jets)):
                        h['MuonElectronJetEvent_deepjet_SingleMuon'].Fill(jets[i].deepjet, genweight)

                    if met_pt > event_cut['metcut'] :

                        h['MuonElectronJethighMetEvent_MuonPt_SingleMuon'].Fill(mu.Pt(), genweight)
                        h['MuonElectronJethighMetEvent_ElectronPt_SingleMuon'].Fill(e.Pt(), genweight)
                        h['MuonElectronJethighMetEvent_JetPt_SingleMuon'].Fill(j.Pt(), genweight)
                        h['MuonElectronJethighMetEvent_dRl_SingleMuon'].Fill(mu.DeltaR(e), genweight)
                        h['MuonElectronJethighMetEvent_dRj_SingleMuon'].Fill(j.DeltaR(e+mu), genweight)
                        h['MuonElectronJethighMetEvent_Mass_SingleMuon'].Fill((mu+e).M(), genweight)
                        for i in range(len(jets)):
                            h['MuonElectronJethighMetEvent_deepjet_SingleMuon'].Fill(jets[i].deepjet, genweight)


    if len(electrons) > 0 and muons[0].charge*electrons[0].charge < 0 :
        
        e = ROOT.TLorentzVector()
        e.SetPtEtaPhiM(electrons[0].pt, electrons[0].eta, electrons[0].phi, electrons[0].mass)

        if ( ( ( mu.Pt() > 8 and e.Pt() > 23 ) or ( mu.Pt() > 23 and e.Pt() > 12 ) ) and isMuonEG == 1 ) : trigger[1] = 1

        if ( mu.Pt() > 50 and isMu == 1 ) or ( mu.Pt() > 27 and isIsoMu == 1 ) : trigger[0] = 1

        if trigger[1] == 1:

            h['MuonElectronEvent_MuonPt_MuonEG'].Fill(mu.Pt(), genweight)
            h['MuonElectronEvent_ElectronPt_MuonEG'].Fill(e.Pt(), genweight)
            h['MuonElectronEvent_MetPt_MuonEG'].Fill(met_pt, genweight)
            h['MuonElectronEvent_dRl_MuonEG'].Fill(mu.DeltaR(e), genweight)
            h['MuonElectronEvent_Mass_MuonEG'].Fill((mu+e).M(), genweight)

            if len(jets) > 0 :

                j = ROOT.TLorentzVector()
                j.SetPtEtaPhiM(jets[0].pt, jets[0].eta, jets[0].phi, jets[0].mass)

                h['MuonElectronJetEvent_MuonPt_MuonEG'].Fill(mu.Pt(), genweight)
                h['MuonElectronJetEvent_ElectronPt_MuonEG'].Fill(e.Pt(), genweight)
                h['MuonElectronJetEvent_JetPt_MuonEG'].Fill(j.Pt(), genweight)
                h['MuonElectronJetEvent_MetPt_MuonEG'].Fill(met_pt, genweight)
                h['MuonElectronJetEvent_dRl_MuonEG'].Fill(mu.DeltaR(e), genweight)
                h['MuonElectronJetEvent_dRj_MuonEG'].Fill(j.DeltaR(e+mu), genweight)
                h['MuonElectronJetEvent_Mass_MuonEG'].Fill((mu+e).M(), genweight)
                for i in range(len(jets)):
                    h['MuonElectronJetEvent_deepjet_MuonEG'].Fill(jets[i].deepjet, genweight)

                if met_pt > event_cut['metcut'] :

                    h['MuonElectronJethighMetEvent_MuonPt_MuonEG'].Fill(mu.Pt(), genweight)
                    h['MuonElectronJethighMetEvent_ElectronPt_MuonEG'].Fill(e.Pt(), genweight)
                    h['MuonElectronJethighMetEvent_JetPt_MuonEG'].Fill(j.Pt(), genweight)
                    h['MuonElectronJethighMetEvent_dRl_MuonEG'].Fill(mu.DeltaR(e), genweight)
                    h['MuonElectronJethighMetEvent_dRj_MuonEG'].Fill(j.DeltaR(e+mu), genweight)
                    h['MuonElectronJethighMetEvent_Mass_MuonEG'].Fill((mu+e).M(), genweight)
                    for i in range(len(jets)):
                        h['MuonElectronJethighMetEvent_deepjet_MuonEG'].Fill(jets[i].deepjet, genweight)


        if ( isData == 0 and ( trigger[0] == 1 or trigger[1] == 1 ) ) or \
           ( isData == 1 and ( isSingleMuonDataset == 1 and trigger[0] == 1 ) ) or \
           ( isData == 1 and ( isMuonEGDataset == 1 and ( trigger[0] == 0 and trigger[1] == 1 ) ) ) :

            if (mu+e).M() > 5.0 :

                h['MuonElectronEvent_MuonPt_Both'].Fill(mu.Pt(), genweight)
                h['MuonElectronEvent_ElectronPt_Both'].Fill(e.Pt(), genweight)
                h['MuonElectronEvent_MetPt_Both'].Fill(met_pt, genweight)
                h['MuonElectronEvent_dRl_Both'].Fill(mu.DeltaR(e), genweight)
                h['MuonElectronEvent_Mass_Both'].Fill((mu+e).M(), genweight)

            if len(jets) > 0 :

                j = ROOT.TLorentzVector()
                j.SetPtEtaPhiM(jets[0].pt, jets[0].eta, jets[0].phi, jets[0].mass)

                if (mu+e).M() > 5.0 and j.Pt() > event_cut['jetPt'] :

                    h['MuonElectronJetEvent_MuonPt_Both'].Fill(mu.Pt(), genweight)
                    h['MuonElectronJetEvent_ElectronPt_Both'].Fill(e.Pt(), genweight)
                    h['MuonElectronJetEvent_JetPt_Both'].Fill(j.Pt(), genweight)
                    h['MuonElectronJetEvent_MetPt_Both'].Fill(met_pt, genweight)
                    h['MuonElectronJetEvent_dRl_Both'].Fill(mu.DeltaR(e), genweight)
                    h['MuonElectronJetEvent_dRj_Both'].Fill(j.DeltaR(e+mu), genweight)
                    h['MuonElectronJetEvent_Mass_Both'].Fill((mu+e).M(), genweight)
                    for i in range(len(jets)):
                        h['MuonElectronJetEvent_deepjet_Both'].Fill(jets[i].deepjet, genweight)

                    if met_pt > event_cut['metcut'] :

                        h['MuonElectronJethighMetEvent_MuonPt_Both'].Fill(mu.Pt(), genweight)
                        h['MuonElectronJethighMetEvent_ElectronPt_Both'].Fill(e.Pt(), genweight)
                        h['MuonElectronJethighMetEvent_JetPt_Both'].Fill(j.Pt(), genweight)
                        h['MuonElectronJethighMetEvent_dRl_Both'].Fill(mu.DeltaR(e), genweight)
                        h['MuonElectronJethighMetEvent_dRj_Both'].Fill(j.DeltaR(e+mu), genweight)
                        for i in range(len(jets)):
                            h['MuonElectronJethighMetEvent_deepjet_Both'].Fill(jets[i].deepjet, genweight)


def Mt(lepton, met):

    cos = np.cos(met.DeltaPhi(lepton))
    Mt = np.sqrt(2*lepton.Pt()*met.Pt()*(1-cos))

    return Mt

def plot_event_hist(region, l1, l2, j, m):
    
    h[region+"_Count"].Fill(0, genweight)

    h[region+"_Mass"].Fill((l1+l2).M(), genweight)
    h[region+"_Lepton1Pt"].Fill(l1.Pt(), genweight)
    h[region+"_Lepton2Pt"].Fill(l2.Pt(), genweight)
    h[region+"_JetPt"].Fill(j.Pt(), genweight)
    h[region+"_MetPt"].Fill(m.Pt(), genweight)
    h[region+"_Nj"].Fill(len(s_j), genweight)
    h[region+"_dRl"].Fill(l1.DeltaR(l2), genweight)
    h[region+"_dRj"].Fill(j.DeltaR(l1+l2), genweight)
    h[region+"_dPhil"].Fill(m.DeltaPhi(l1), m.DeltaPhi(l2), genweight)
    h[region+"_dPhi"].Fill(m.DeltaPhi(l1+l2), m.DeltaPhi(j), genweight)

    h[region+"_Mtl1"].Fill(Mt(l1,m), genweight)
    h[region+"_Mtl2"].Fill(Mt(l2,m), genweight)
    h[region+"_Mtl"].Fill(Mt((l1+l2),m), genweight)
    
    h[region+"_cosMl1"].Fill(np.cos(m.DeltaPhi(l1)), genweight)
    h[region+"_cosMl2"].Fill(np.cos(m.DeltaPhi(l2)), genweight)

    h[region+"_cosl"].Fill(np.cos(l1.DeltaPhi(l2)), genweight)

    if len(s_b) > 0 :
        for i in s_b:
            h[region+"_bJetPt"].Fill(i.pt, genweight)
    

inputFileNames=open(inputFileList, 'r')
for inputFileName in inputFileNames:
    inputFileName=inputFileName.replace("\n","")
    print(inputFileName.replace("\n",""))

    fchain.Add(inputFileName)
    chain2.Add(inputFileName)
    if isData == 0:
        chain3.Add(inputFileName)
        chain4.Add(inputFileName)

fchain.AddFriend(chain2)
if isData == 0:
    fchain.AddFriend(chain3)
    fchain.AddFriend(chain4)

jets = ROOT.JetInfoDS()
muons = ROOT.MuonInfoDS()
electrons = ROOT.ElectronInfoDS()
tausUnCleaned = ROOT.TauInfoDS()
tausECleaned = ROOT.TauInfoDS()
tausMCleaned = ROOT.TauInfoDS()
tausBoosted = ROOT.TauInfoDS()

fchain.SetBranchAddress("Jets", ROOT.AddressOf(jets))
fchain.SetBranchAddress("Muons", ROOT.AddressOf(muons))
fchain.SetBranchAddress("Electrons", ROOT.AddressOf(electrons))
fchain.SetBranchAddress("TausUnCleaned", ROOT.AddressOf(tausUnCleaned))
fchain.SetBranchAddress("TausECleaned", ROOT.AddressOf(tausECleaned))
fchain.SetBranchAddress("TausMCleaned", ROOT.AddressOf(tausMCleaned))
fchain.SetBranchAddress("TausBoosted", ROOT.AddressOf(tausBoosted))

if isData == 0:
    genParticle = ROOT.GenParticleInfoDS()
    fchain.SetBranchAddress("GenParticleInfo", ROOT.AddressOf(genParticle))

for iev in range(fchain.GetEntries()): # Be careful!!!                                                               

   fchain.GetEntry(iev)

   mets = fchain.GetBranch("Mets")
   met_pt = mets.GetLeaf('pt').GetValue()
   met_phi = mets.GetLeaf('phi').GetValue()

   h['hEvents'].Fill(0.5, 1)

   h['hMetPt'].Fill(met_pt)

   if isData == 1: genweight = 1

   if isData == 0:
       weight = fchain.GetBranch("lumiInfo")
       genweight = weight.GetLeaf('weight').GetValue()

       s_genmu = []
       
       if genParticle.size()>0:
           for i in range(genParticle.size()):
               gen = genParticle.at(i)
               if abs(gen.pdgid) == 13 : s_genmu += [gen]

   h['hEvents'].Fill(1.5, genweight)

   isSingleJet = fchain.GetLeaf('isSingleJet').GetValue()
   isHT = fchain.GetLeaf('isHT').GetValue()
   isHTMHT = fchain.GetLeaf('isHTMHT').GetValue()
   isMu = fchain.GetLeaf('isMu').GetValue()
   isIsoMu = fchain.GetLeaf('isIsoMu').GetValue()
   isIsoMuTau = fchain.GetLeaf('isIsoMuTau').GetValue()
   isIsoEle = fchain.GetLeaf('isIsoEle').GetValue()
   isEleTau = fchain.GetLeaf('isEleTau').GetValue()
   isMuonEG = fchain.GetLeaf('isMuonEG').GetValue()

   h["isSingleJet"].Fill(isSingleJet, genweight)
   h["isHT"].Fill(isHT, genweight)

   s_j = []
   s_b = []
   s_e = []
   s_lowE = []
   s_mu = []
   s_isomu = []
   s_isoe = []

   s_nonIsoMu = []
   s_nonIsoE = []

   unclean = []
   eclean = []
   mclean = []
   lowEclean = []

   boosted = []
   mclean_altered = []

   if jets.size()>0:
      for i in range(jets.size()):
         jet = jets.at(i)
         if abs(jet.eta) < 2.5 :
             if jet.id >= 2:
                 h['hJetPt'].Fill(jet.pt, genweight)
                 h['hBtag'].Fill(jet.deepjet, genweight)
                 s_j+=[jet]
                 if jet.deepjet > 0.7476:
                     h['hBJetPt'].Fill(jet.pt, genweight)
                     s_b+=[jet]

   if muons.size()>0:
       for i in range(muons.size()):
           muon = muons.at(i)
           if abs(muon.eta) < 2.4 and muon.pt > 8.0:
               if muon.id >= 1:
                   h['hMuonPt'].Fill(muon.pt, genweight) 
                   s_mu+=[muon]
                   if muon.iso < 0.25:
                       h['hIsoMuonPt'].Fill(muon.pt, genweight)
                       s_isomu+=[muon]
                   if muon.iso > 0.25:
                       h['hNonIsoMuonPt'].Fill(muon.pt, genweight)
                       s_nonIsoMu+=[muon]

   if electrons.size()>0:
      for i in range(electrons.size()):
         electron = electrons.at(i)
         if abs(electron.eta) < 2.5 and electron.pt > 12.0 :
             if electron.id >= 1 :
                 h['hElectronPt'].Fill(electron.pt, genweight)
                 s_e+=[electron]
                 if electron.iso >= 1:
                     h['hIsoElectronPt'].Fill(electron.pt, genweight)
                     s_isoe+=[electron]
                 if electron.iso == 0:
                     h['hNonIsoElectronPt'].Fill(electron.pt, genweight)
                     s_nonIsoE+=[electron]


   s_j.sort(key=lambda x: x.pt, reverse=True)
   s_e.sort(key=lambda x: x.pt, reverse=True)
   s_isoe.sort(key=lambda x: x.pt, reverse=True)
   s_mu.sort(key=lambda x: x.pt, reverse=True)
   s_isomu.sort(key=lambda x: x.pt, reverse=True)

   s_nonIsoMu.sort(key=lambda x: x.pt, reverse=True)
   s_nonIsoE.sort(key=lambda x: x.pt, reverse=True)

   isEMu = 0
   isMuMu = 0



   if len(s_isomu) == 1 and len(s_isoe) == 1 and len(s_j) > 1 and s_isomu[0].charge*s_isoe[0].charge < 0 :
       kin_selection(s_isomu, s_isoe, s_j)

   if len(s_isomu) > 0 :
       check_cut_flow(s_isomu, s_isoe, s_j)
   
   if len(s_isomu) > 1 and len(s_j) > 0 and s_isomu[0].charge*s_isomu[1].charge < 0 : 
       MuMu_Channel(s_isomu, s_j, met_pt, met_phi, plot=True)

   if len(s_isomu) > 0 and len(s_isoe) > 0 and len(s_j) > 0 and s_isoe[0].charge*s_isomu[0].charge < 0 :
       EMu_Channel(s_isoe, s_isomu, s_j, met_pt, met_phi,)

   if len(s_isomu) > 1 and len(s_j) > 0 and s_isomu[0].charge*s_isomu[1].charge < 0 and len(s_b) == 0 :
       if MuMu_Channel(s_isomu, s_j, met_pt, met_phi) == 1 : continue

   if len(s_isomu) > 0 and len(s_isoe) > 0 and len(s_j) > 0 and s_isoe[0].charge*s_isomu[0].charge < 0 :
       EMu_Channel(s_isoe, s_isomu, s_j, met_pt, met_phi, "hEMu_checkMuMu")

   if isNonIso == True :
       if len(s_nonIsoMu) > 0 and len(s_isoe) > 0 and len(s_j) > 0 and s_nonIsoMu[0].charge*s_isoe[0].charge < 0 :
           EMu_Channel(s_isoe, s_nonIsoMu, s_j, met_pt, met_phi, "hEMu_nonIsoMu")

       if len(s_isomu) > 0 and len(s_nonIsoE) > 0 and len(s_j) > 0 and s_isomu[0].charge*s_nonIsoE[0].charge < 0 :
           EMu_Channel(s_nonIsoE, s_isomu, s_j, met_pt, met_phi, "hEMu_nonIsoE")

       if len(s_mu) > 0 and len(s_e) > 0 and len(s_j) > 0 and s_mu[0].charge*s_e[0].charge < 0 :
           EMu_Channel(s_e, s_mu, s_j, met_pt, met_phi, "hEMu_nonIso")


out.cd()

for key in h.keys():
    h[key].Write()

out.Close()

print("--- %s seconds ---" % (time.time() - start_time))
