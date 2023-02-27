import ROOT, sys, os
import numpy as np
import time

start_time = time.time()

opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]

plotSemiLeptonic = 0

plotSingleMuonHTTrig = 0
plotJetHT = 0
plotSingleMuon = 0
plotHTTrig = 0

plot2DforTau = 0
isAltered = 0

print("isAltered = ", isAltered)
print("plot2DforTau = ", plot2DforTau)

if "-b" in opts : 
    isData = 0
    plotJetHT = 0
    plotSingleMuon = 1
    plotHTTrig = 1
    if plotSingleMuon == 1 and plotHTTrig == 1 :
        plotSingleMuonHTTrig = 1

isJetHTDataset = 0
isSingleMuonDataset = 0

if "-dj" in opts :
    isData = 1
    isJetHTDataset = 1
    isSingleMuonDataset = 0

if "-ds" in opts :
    isData = 1
    isJetHTDataset = 0
    isSingleMuonDataset = 1

print("isData = ", isData)
print("isJetHTDataset = ", isJetHTDataset)
print("isSingleMuonDataset = ", isSingleMuonDataset)

if isJetHTDataset == 1 :
    plotHTTrig = 1
    plotSingleMuon = 0
    
if isSingleMuonDataset == 1 :
    plotSingleMuon = 1

ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/JetInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/MuonInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/ElectronInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/TauInfoDS.h"')

inputFileListName=sys.argv[1]
inputFileList=inputFileListName

print("Parameters:")
print("plotSemiLeptonic = ",plotSemiLeptonic)
print("plotJetHT = ",plotJetHT)
print("plotHTTrig = ",plotHTTrig)
print("plotSingleMuon = ",plotSingleMuon)
print("plotSingleMuonHTTrig = ",plotSingleMuonHTTrig)

if len(sys.argv)>2 and not sys.argv[2].startswith("-"):
    outputFileDir=sys.argv[2]
else:
    outputFileDir = "./output/"


if plotSingleMuonHTTrig == 1 and isData== 0 :
    if plotSemiLeptonic == 1:
        outputFileName = outputFileDir+"h_debugMuTau_HighHT_FullyLeptonic_Inclusive_MCOnly_"+inputFileListName.split("/")[-1].replace(".txt",".root")
        if isAltered == 1 :
            outputFileName = outputFileDir+"h_debugMuTau_HighHT_FullyLeptonic_Inclusive_Altered_MCOnly_"+inputFileListName.split("/")[-1].replace(".txt",".root")
    elif plotSemiLeptonic == 0:
        outputFileName = outputFileDir+"h_debugMuTau_HighHT_Inclusive_MCOnly_"+inputFileListName.split("/")[-1].replace(".txt",".root")
        if isAltered == 1 :
            outputFileName = outputFileDir+"h_debugMuTau_HighHT_Inclusive_Altered_MCOnly_"+inputFileListName.split("/")[-1].replace(".txt",".root")


out=ROOT.TFile.Open(outputFileName,'recreate')
print(outputFileName)

fchain = ROOT.TChain('tcpNtuples/analysisTree')
chain2 = ROOT.TChain('tcpTrigNtuples/triggerTree')
if isData == 0:
    chain3 = ROOT.TChain('lumiSummary/lumiTree')
    chain4 = ROOT.TChain('tcpGenNtuples/genTree')

pi = np.pi


def book1DHist(region):

    h[region+"_Mass"] = ROOT.TH1F (region+"_Mass", region+"_Mass ; M_{vis.} (GeV) ; Events ", 100, 0, 100)
    h[region+"_Lepton1Pt"] = ROOT.TH1F (region+"_Lepton1Pt", region+"_Lepton1Pt ; P_{T} (GeV) ; Events ", 500, 0, 500)
    h[region+"_Lepton2Pt"] = ROOT.TH1F (region+"_Lepton2Pt", region+"_Lepton2Pt ; P_{T} (GeV) ; Events ", 500, 0, 500)
    h[region+"_JetPt"] = ROOT.TH1F (region+"_JetPt", region+"_JetPt ; JetP_{T} (GeV) ; Events ", 2000, 0, 2000)
    h[region+"_Mt"] = ROOT.TH1F (region+"_Mt", region+"_Mt ; M_{T} (GeV) ; Events ", 150, 0, 150)
    h[region+"_MetPt"] = ROOT.TH1F (region+"_MetPt", region+"_MetPt ; MET (GeV) ; Events ", 500, 0, 500)
    h[region+"_Nj"] = ROOT.TH1F (region+"_Nj", region+"_Nj ; N_{j} ; Events ", 10, 0, 10)


h = {}

h['ChOverlap'] = ROOT.TH1F("ChOverlap", "Channel Overlap ; ChOverlap; N", 12, 1, 13)

h['isSingleJet'] = ROOT.TH1F("isSingleJet", "isSingleJet ; isSingleJet ; N", 4,-1.5,2.5)
h['isHT'] = ROOT.TH1F("isHT", "isHT ; isHT ; N", 4,-1.5,2.5)

h['hEvents'] = ROOT.TH1F ("NEvents", "Number of Events; ;N", 2, 0, 2)
h['hTauEvents'] = ROOT.TH1F ("hTauEvents", "Number of each decay modes;;a.u.",10,0,10)

h['hMuMu_Events'] = ROOT.TH1F ("hMuMu_Events", "hMuMu_Events;;N", 6, 1, 7)
h['hEMu_Events'] = ROOT.TH1F ("hEMu_Events", "hEMu_Events;;N", 6, 1, 7)
h['hMuTau_Events'] = ROOT.TH1F ("hMuTau_Events", "hMuTau_Events;;N", 6, 1, 7)
h['hMuTau_2DPlot_Events'] = ROOT.TH1F ("hMuTau_2DPlot_Events", "hMuTau_2DPlot_Events;;N", 6, 1, 7)
h['hMuTau_SR_Events'] = ROOT.TH1F ("hMuTau_SR_Events", "hMuTau_SR_Events;;N", 6, 1, 7)
h['hMuTau_SR_2DPlot_Events'] = ROOT.TH1F ("hMuTau_SR_2DPlot_Events", "hMuTau_SR_2DPlot_Events;;N", 6, 1, 7)

h['hMuTau_DR1_Events'] = ROOT.TH1F ("hMuTau_DR1_Events", "hMuTau_DR1_Events;;N", 4, 0, 4)
h['hMuTau_DR2_Events'] = ROOT.TH1F ("hMuTau_DR2_Events", "hMuTau_DR2_Events;;N", 4, 0, 4)
h['hMuTau_DR3_Events'] = ROOT.TH1F ("hMuTau_DR3_Events", "hMuTau_DR3_Events;;N", 4, 0, 4)
h['hMuTau_DR4_Events'] = ROOT.TH1F ("hMuTau_DR4_Events", "hMuTau_DR4_Events;;N", 4, 0, 4)
h['hMuTau_DR5_Events'] = ROOT.TH1F ("hMuTau_DR5_Events", "hMuTau_DR5_Events;;N", 4, 0, 4)
h['hMuTau_DR6_Events'] = ROOT.TH1F ("hMuTau_DR6_Events", "hMuTau_DR6_Events;;N", 4, 0, 4)
h['hMuTau_DR7_Events'] = ROOT.TH1F ("hMuTau_DR7_Events", "hMuTau_DR7_Events;;N", 4, 0, 4)

h['hMuMu_Trigger'] = ROOT.TH1F ("hMuMu_Trigger", "hMuMu_Trigger;;N", 7, 1, 8)
h['hEMu_Trigger'] = ROOT.TH1F ("hEMu_Trigger", "hEMu_Trigger;;N", 7, 1, 8)
h['hMuTau_Trigger'] = ROOT.TH1F ("hMuTau_Trigger", "hMuTau_Trigger;;N", 7, 1, 8)
h['hMuTau_SR_Trigger'] = ROOT.TH1F ("hMuTau_SR_Trigger", "hMuTau_SR_Trigger;;N", 7, 1, 8)

h['hMuMu_Inclusive_Trigger'] = ROOT.TH1F ("hMuMu_Inclusive_Trigger", "hMuMu_Inclusive_Trigger;;N", 5, 1, 6)
h['hEMu_Inclusive_Trigger'] = ROOT.TH1F ("hEMu_Inclusive_Trigger", "hEMu_Inclusive_Trigger;;N", 5, 1, 6)
h['hMuTau_Inclusive_Trigger'] = ROOT.TH1F ("hMuTau_Inclusive_Trigger", "hMuTau_Inclusive_Trigger;;N", 5, 1, 6)
h['hMuTau_SR_Inclusive_Trigger'] = ROOT.TH1F ("hMuTau_SR_Inclusive_Trigger", "hMuTau_SR_Inclusive_Trigger;;N", 7, 1, 8)

h['hMuMu_SR_dRcut_highMET_dPhicut'] = ROOT.TH1F ("hMuMu_SR_dRcut_highMET_dPhicut", "hMuMu_SR_dRcut_highMET_dPhicut; ;N", 100, 0, 100)
h['hEMu_SR_dRcut_highMET_dPhicut'] = ROOT.TH1F ("hEMu_SR_dRcut_highMET_dPhicut", "hEMu_SR_dRcut_highMET_dPhicut; ;N", 100, 0, 100)

if plotSemiLeptonic == 1:
    h['hMuMu_Baseline_Nbj'] = ROOT.TH1F ("hMuMu_Baseline_Nbj", "hMuMu_Baseline_Nbj;;N", 5, 0, 5)
    h['hMuMu_highMET_Nbj'] = ROOT.TH1F ("hMuMu_highMET_Nbj", "hMuMu_highMET_Nbj;;N", 5, 0, 5)
    h['hMuMu_dRcut_Nbj'] = ROOT.TH1F ("hMuMu_dRcut_Nbj", "hMuMu_dRcut_Nbj;;N", 5, 0, 5)

h['hMuTau_OS_Nbj'] = ROOT.TH1F ("hMuTau_OS_Nbj", "hMuTau_OS_Nbj;;N", 5, 0, 5)
h['hMuTau_SS_Nbj'] = ROOT.TH1F ("hMuTau_SS_Nbj", "hMuTau_SS_Nbj;;N", 5, 0, 5)

if plotSemiLeptonic == 1:
    h['hMuMu_lowMET'] = ROOT.TH1F ("hMuMu_lowMET", "hMuMu_lowMET; ;N", 100, 0, 100)
    h['hMuMu_lowMET_dRl'] = ROOT.TH1F ("hMuMu_lowMET_dRl", "hMuMu_lowMET_dRl;;N", 100, 0, 5)
    h['hMuMu_lowMET_dRj'] = ROOT.TH1F ("hMuMu_lowMET_dRj", "hMuMu_lowMET_dRj;;N", 100, 0, 5)
    h['hMuMu_lowMET_MetPt'] = ROOT.TH1F ("hMuMu_lowMET_MetPt", "hMuMu_lowMET_MetPt;;N", 100, 0, 100)
    h['hMuMu_lowMET_Nj'] = ROOT.TH1F ("hMuMu_lowMET_Nj", "hMuMu_lowMET_Nj;;N", 10, 0, 10)
    h['hMuMu_lowMET_Muon1Pt'] = ROOT.TH1F ("hMuMu_lowMET_Muon1Pt", "hMuMu_lowMET_Muon1Pt;;N", 500, 0, 500)
    h['hMuMu_lowMET_Muon2Pt'] = ROOT.TH1F ("hMuMu_lowMET_Muon2Pt", "hMuMu_lowMET_Muon2Pt;;N", 500, 0, 500)
    h['hMuMu_lowMET_JetPt'] = ROOT.TH1F ("hMuMu_lowMET_JetPt", "hMuMu_lowMET_JetPt;;N", 2000, 0, 2000)
    h['hMuMu_lowMET_dPhil'] = ROOT.TH1F ("hMuMu_lowMET_dPhil", "hMuMu_lowMET_dPhil;;N", 100, -pi, pi)
    h['hMuMu_lowMET_dPhij'] = ROOT.TH1F ("hMuMu_lowMET_dPhij", "hMuMu_lowMET_dPhij;;N", 100, -pi, pi)

    h['hMuMu_dRcut'] = ROOT.TH1F ("hMuMu_dRcut", "hMuMu_dRcut; ;N", 100, 0, 100)
    h['hMuMu_dRcut_dRl'] = ROOT.TH1F ("hMuMu_dRcut_dRl", "hMuMu_dRcut_dRl;;N", 20, 0, 0.5)
    h['hMuMu_dRcut_dRj'] = ROOT.TH1F ("hMuMu_dRcut_dRj", "hMuMu_dRcut_dRj;;N", 100, 0, 5)
    h['hMuMu_dRcut_MetPt'] = ROOT.TH1F ("hMuMu_dRcut_MetPt", "hMuMu_dRcut_MetPt;;N", 100, 0, 100)
    h['hMuMu_dRcut_Nj'] = ROOT.TH1F ("hMuMu_dRcut_Nj", "hMuMu_dRcut_Nj;;N", 5, 0, 5)
    h['hMuMu_lowMET_dRcut_Nbj'] = ROOT.TH1F ("hMuMu_lowMET_dRcut_Nbj", "hMuMu_lowMET_dRcut_Nbj;;N", 10, 0, 10)
    h['hMuMu_dRcut_Muon1Pt'] = ROOT.TH1F ("hMuMu_dRcut_Muon1Pt", "hMuMu_dRcut_Muon1Pt;;N", 500, 0, 500)
    h['hMuMu_dRcut_Muon2Pt'] = ROOT.TH1F ("hMuMu_dRcut_Muon2Pt", "hMuMu_dRcut_Muon2Pt;;N", 500, 0, 500)
    h['hMuMu_dRcut_JetPt'] = ROOT.TH1F ("hMuMu_dRcut_JetPt", "hMuMu_dRcut_JetPt;;N", 2000, 0, 2000)
    h['hMuMu_dRcut_dPhil'] = ROOT.TH1F ("hMuMu_dRcut_dPhil", "hMuMu_dRcut_dPhil;;N", 100, -pi, pi)
    h['hMuMu_dRcut_dPhij'] = ROOT.TH1F ("hMuMu_dRcut_dPhij", "hMuMu_dRcut_dPhij;;N", 100, -pi, pi)

    h['hMuMu_dPhicut'] = ROOT.TH1F ("hMuMu_dPhicut", "hMuMu_dPhicut; ;N", 100, 0, 100)
    h['hMuMu_dPhicut_dRl'] = ROOT.TH1F ("hMuMu_dPhicut_dRl", "hMuMu_dPhicut_dRl;;N", 100, 0, 5)
    h['hMuMu_dPhicut_dRj'] = ROOT.TH1F ("hMuMu_dPhicut_dRj", "hMuMu_dPhicut_dRj;;N", 100, 0, 5)
    h['hMuMu_dPhicut_MetPt'] = ROOT.TH1F ("hMuMu_dPhicut_MetPt", "hMuMu_dPhicut_MetPt;;N", 100, 0, 100)
    h['hMuMu_dPhicut_Nj'] = ROOT.TH1F ("hMuMu_dPhicut_Nj", "hMuMu_dPhicut_Nj;;N", 10, 0, 10)
    h['hMuMu_dPhicut_Nbj'] = ROOT.TH1F ("hMuMu_dPhicut_Nbj", "hMuMu_dPhicut_Nbj;;N", 5, 0, 5)
    h['hMuMu_dPhicut_Muon1Pt'] = ROOT.TH1F ("hMuMu_dPhicut_Muon1Pt", "hMuMu_dPhicut_Muon1Pt;;N", 500, 0, 500)
    h['hMuMu_dPhicut_Muon2Pt'] = ROOT.TH1F ("hMuMu_dPhicut_Muon2Pt", "hMuMu_dPhicut_Muon2Pt;;N", 500, 0, 500)
    h['hMuMu_dPhicut_JetPt'] = ROOT.TH1F ("hMuMu_dPhicut_JetPt", "hMuMu_dPhicut_JetPt;;N", 2000, 0, 2000)
    h['hMuMu_dPhicut_dPhil'] = ROOT.TH1F ("hMuMu_dPhicut_dPhil", "hMuMu_dPhicut_dPhil;;N", 100, -pi, pi)
    h['hMuMu_dPhicut_dPhij'] = ROOT.TH1F ("hMuMu_dPhicut_dPhij", "hMuMu_dPhicut_dPhij;;N", 100, -pi, pi)

#############-------------------------------EMu

if plotSemiLeptonic == 1:
    h['hEMu_lowMET'] = ROOT.TH1F ("hEMu_lowMET", "hEMu_lowMET; ;N", 100, 0, 100)
    h['hEMu_lowMET_dRl'] = ROOT.TH1F ("hEMu_lowMET_dRl", "hEMu_lowMET_dRl;;N", 100, 0, 5)
    h['hEMu_lowMET_dRj'] = ROOT.TH1F ("hEMu_lowMET_dRj", "hEMu_lowMET_dRj;;N", 100, 0, 5)
    h['hEMu_lowMET_MetPt'] = ROOT.TH1F ("hEMu_lowMET_MetPt", "hEMu_lowMET_MetPt;;N", 100, 0, 100)
    h['hEMu_lowMET_Nj'] = ROOT.TH1F ("hEMu_lowMET_Nj", "hEMu_lowMET_Nj;;N", 10, 0, 10)
    h['hEMu_lowMET_MuonPt'] = ROOT.TH1F ("hEMu_lowMET_MuonPt", "hEMu_lowMET_MuonPt;;N", 500, 0, 500)
    h['hEMu_lowMET_ElectronPt'] = ROOT.TH1F ("hEMu_lowMET_ElectronPt", "hEMu_lowMET_ElectronPt;;N", 500, 0, 500)
    h['hEMu_lowMET_JetPt'] = ROOT.TH1F ("hEMu_lowMET_JetPt", "hEMu_lowMET_JetPt;;N", 2000, 0, 2000)
    h['hEMu_lowMET_dPhil'] = ROOT.TH1F ("hEMu_lowMET_dPhil", "hEMu_lowMET_dPhil;;N", 100, -pi, pi)
    h['hEMu_lowMET_dPhij'] = ROOT.TH1F ("hEMu_lowMET_dPhij", "hEMu_lowMET_dPhij;;N", 100, -pi, pi)

    h['hEMu_dRcut'] = ROOT.TH1F ("hEMu_dRcut", "hEMu_dRcut; ;N", 100, 0, 100)
    h['hEMu_dRcut_dRl'] = ROOT.TH1F ("hEMu_dRcut_dRl", "hEMu_dRcut_dRl;;N", 20, 0, 0.5)
    h['hEMu_dRcut_dRj'] = ROOT.TH1F ("hEMu_dRcut_dRj", "hEMu_dRcut_dRj;;N", 100, 0, 5)
    h['hEMu_dRcut_MetPt'] = ROOT.TH1F ("hEMu_dRcut_MetPt", "hEMu_dRcut_MetPt;;N", 100, 0, 100)
    h['hEMu_dRcut_Nj'] = ROOT.TH1F ("hEMu_dRcut_Nj", "hEMu_dRcut_Nj;;N", 10, 0, 10)
    h['hEMu_lowMET_dRcut_Nbj'] = ROOT.TH1F ("hEMu_lowMET_dRcut_Nbj", "hEMu_lowMET_dRcut_Nbj;;N", 5, 0, 5)
    h['hEMu_dRcut_MuonPt'] = ROOT.TH1F ("hEMu_dRcut_MuonPt", "hEMu_dRcut_MuonPt;;N", 500, 0, 500)
    h['hEMu_dRcut_ElectronPt'] = ROOT.TH1F ("hEMu_dRcut_ElectronPt", "hEMu_dRcut_ElectronPt;;N", 500, 0, 500)
    h['hEMu_dRcut_JetPt'] = ROOT.TH1F ("hEMu_dRcut_JetPt", "hEMu_dRcut_JetPt;;N", 2000, 0, 2000)
    h['hEMu_dRcut_dPhil'] = ROOT.TH1F ("hEMu_dRcut_dPhil", "hEMu_dRcut_dPhil;;N", 100, -pi, pi)
    h['hEMu_dRcut_dPhij'] = ROOT.TH1F ("hEMu_dRcut_dPhij", "hEMu_dRcut_dPhij;;N", 100, -pi, pi)

#############-------------------------------MuTau

if plot2DforTau == 0:

    h['hMuTau_lowMET_dRl'] = ROOT.TH1F ("hMuTau_lowMET_dRl", "hMuTau_lowMET_dRl;;N", 100, 0, 5)
    h['hMuTau_lowMET_dRj'] = ROOT.TH1F ("hMuTau_lowMET_dRj", "hMuTau_lowMET_dRj;;N", 100, 0, 5)

    h['hMuTau_lowMET_highMt_dRl'] = ROOT.TH1F ("hMuTau_lowMET_highMt_dRl", "hMuTau_lowMET_highMt_dRl;;N", 100, 0, 5)
    h['hMuTau_lowMET_highMt_dRj'] = ROOT.TH1F ("hMuTau_lowMET_highMt_dRj", "hMuTau_lowMET_highMt_dRj;;N", 100, 0, 5)

    h['hMuTau_lowMET_lowMt_dRl'] = ROOT.TH1F ("hMuTau_lowMET_lowMt_dRl", "hMuTau_lowMET_lowMt_dRl;;N", 100, 0, 5)
    h['hMuTau_lowMET_lowMt_dRj'] = ROOT.TH1F ("hMuTau_lowMET_lowMt_dRj", "hMuTau_lowMET_lowMt_dRj;;N", 100, 0, 5)

    h['hMuTau_lowMET_dRcut_lowMt_TauPt0'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_lowMt_TauPt0", "hMuTau_lowMET_dRcut_lowMt_TauPt0;;N", 500, 0, 500)
    h['hMuTau_lowMET_dRcut_lowMt_TauPt1'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_lowMt_TauPt1", "hMuTau_lowMET_dRcut_lowMt_TauPt1;;N", 500, 0, 500)
    h['hMuTau_lowMET_dRcut_lowMt_TauPt10'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_lowMt_TauPt10", "hMuTau_lowMET_dRcut_lowMt_TauPt10;;N", 500, 0, 500)

    h['hMuTau_lowMET_dRcut_highMt_TauPt0'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_highMt_TauPt0", "hMuTau_lowMET_dRcut_highMt_TauPt0;;N", 500, 0, 500)
    h['hMuTau_lowMET_dRcut_highMt_TauPt1'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_highMt_TauPt1", "hMuTau_lowMET_dRcut_highMt_TauPt1;;N", 500, 0, 500)
    h['hMuTau_lowMET_dRcut_highMt_TauPt10'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_highMt_TauPt10", "hMuTau_lowMET_dRcut_highMt_TauPt10;;N", 500, 0, 500)

    h['hMuTau_highMt_dRl'] = ROOT.TH1F ("hMuTau_highMt_dRl", "hMuTau_highMt_dRl;;N", 100, 0, 5)
    h['hMuTau_highMt_dRj'] = ROOT.TH1F ("hMuTau_highMt_dRj", "hMuTau_highMt_dRj;;N", 100, 0, 5)
    h['hMuTau_highMt_highMET_dRl'] = ROOT.TH1F ("hMuTau_highMt_highMET_dRl", "hMuTau_highMt_highMET_dRl;;N", 100, 0, 5)
    h['hMuTau_highMt_highMET_dRj'] = ROOT.TH1F ("hMuTau_highMt_highMET_dRj", "hMuTau_highMt_highMET_dRj;;N", 100, 0, 5)

    h['hMuTau_highMt_dRcut_highMET_TauPt0'] = ROOT.TH1F ("hMuTau_highMt_dRcut_highMET_TauPt0", "hMuTau_highMt_dRcut_highMET_TauPt0;;N", 500, 0, 500)
    h['hMuTau_highMt_dRcut_highMET_TauPt1'] = ROOT.TH1F ("hMuTau_highMt_dRcut_highMET_TauPt1", "hMuTau_highMt_dRcut_highMET_TauPt1;;N", 500, 0, 500)
    h['hMuTau_highMt_dRcut_highMET_TauPt10'] = ROOT.TH1F ("hMuTau_highMt_dRcut_highMET_TauPt10", "hMuTau_highMt_dRcut_highMET_TauPt10;;N", 500, 0, 500)

    h['hMuTau_OS_TauPt'] = ROOT.TH1F ("hMuTau_OS_TauPt", "hMuTau_OS_TauPt;;N", 500, 0, 500)
    h['hMuTau_OS'] = ROOT.TH1F ("hMuTau_OS", "hMuTau_OS;;N", 100, 0, 100)
    h['hMuTau_OS_Nj'] = ROOT.TH1F ("hMuTau_OS_Nj", "hMuTau_OS_Nj;;N", 10, 0, 10)

    h['hMuTau_SS_dRl'] = ROOT.TH1F ("hMuTau_SS_dRl", "hMuTau_SS_dRl;;N", 100, 0, 5)
    h['hMuTau_SS_dRj'] = ROOT.TH1F ("hMuTau_SS_dRj", "hMuTau_SS_dRj;;N", 100, 0, 5)

    h['hMuTau_SS_lowMET_dRl'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRl", "hMuTau_SS_lowMET_dRl;;N", 100, 0, 5)
    h['hMuTau_SS_lowMET_dRj'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRj", "hMuTau_SS_lowMET_dRj;;N", 100, 0, 5)

    h['hMuTau_SS_lowMET_highMt_dRl'] = ROOT.TH1F ("hMuTau_SS_lowMET_highMt_dRl", "hMuTau_SS_lowMET_highMt_dRl;;N", 100, 0, 5)
    h['hMuTau_SS_lowMET_highMt_dRj'] = ROOT.TH1F ("hMuTau_SS_lowMET_highMt_dRj", "hMuTau_SS_lowMET_highMt_dRj;;N", 100, 0, 5)

    h['hMuTau_SS_lowMET_lowMt_dRl'] = ROOT.TH1F ("hMuTau_SS_lowMET_lowMt_dRl", "hMuTau_SS_lowMET_lowMt_dRl;;N", 100, 0, 5)
    h['hMuTau_SS_lowMET_lowMt_dRj'] = ROOT.TH1F ("hMuTau_SS_lowMET_lowMt_dRj", "hMuTau_SS_lowMET_lowMt_dRj;;N", 100, 0, 5)

    h['hMuTau_SS_lowMET_dRcut_highMt_TauPt0'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_highMt_TauPt0", "hMuTau_SS_lowMET_dRcut_highMt_TauPt0;;N", 500, 0, 500)
    h['hMuTau_SS_lowMET_dRcut_highMt_TauPt1'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_highMt_TauPt1", "hMuTau_SS_lowMET_dRcut_highMt_TauPt1;;N", 500, 0, 500)
    h['hMuTau_SS_lowMET_dRcut_highMt_TauPt10'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_highMt_TauPt10", "hMuTau_SS_lowMET_dRcut_highMt_TauPt10;;N", 500, 0, 500)

    h['hMuTau_SS_lowMET_dRcut_lowMt_TauPt0'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_lowMt_TauPt0", "hMuTau_SS_lowMET_dRcut_lowMt_TauPt0;;N", 500, 0, 500)
    h['hMuTau_SS_lowMET_dRcut_lowMt_TauPt1'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_lowMt_TauPt1", "hMuTau_SS_lowMET_dRcut_lowMt_TauPt1;;N", 500, 0, 500)
    h['hMuTau_SS_lowMET_dRcut_lowMt_TauPt10'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_lowMt_TauPt10", "hMuTau_SS_lowMET_dRcut_lowMt_TauPt10;;N", 500, 0, 500)

    h['hMuTau_SS_dRcut_highMET_highMt_TauPt0'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_highMt_TauPt0", "hMuTau_SS_dRcut_highMET_highMt_TauPt0;;N", 500, 0, 500)
    h['hMuTau_SS_dRcut_highMET_highMt_TauPt1'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_highMt_TauPt1", "hMuTau_SS_dRcut_highMET_highMt_TauPt1;;N", 500, 0, 500)
    h['hMuTau_SS_dRcut_highMET_highMt_TauPt10'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_highMt_TauPt10", "hMuTau_SS_dRcut_highMET_highMt_TauPt10;;N", 500, 0, 500)

    h['hMuTau_SS_dRcut_highMET_lowMt_TauPt0'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_lowMt_TauPt0", "hMuTau_SS_dRcut_highMET_lowMt_TauPt0;;N", 500, 0, 500)
    h['hMuTau_SS_dRcut_highMET_lowMt_TauPt1'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_lowMt_TauPt1", "hMuTau_SS_dRcut_highMET_lowMt_TauPt1;;N", 500, 0, 500)
    h['hMuTau_SS_dRcut_highMET_lowMt_TauPt10'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_lowMt_TauPt10", "hMuTau_SS_dRcut_highMET_lowMt_TauPt10;;N", 500, 0, 500)

    if isData == 0 or isAltered == 1:
        h['hMuTau_SR_dRl'] = ROOT.TH1F ("hMuTau_SR_dRl", "hMuTau_SR_dRl;;N", 100, 0, 5)
        h['hMuTau_SR_dRj'] = ROOT.TH1F ("hMuTau_SR_dRj", "hMuTau_SR_dRj;;N", 100, 0, 5)

        h['hMuTau_SR_dRcut_dRl'] = ROOT.TH1F ("hMuTau_SR_dRcut_dRl", "hMuTau_SR_dRcut_dRl;;N", 100, 0, 5)
        h['hMuTau_SR_dRcut_dRj'] = ROOT.TH1F ("hMuTau_SR_dRcut_dRj", "hMuTau_SR_dRcut_dRj;;N", 100, 0, 5)

        h['hMuTau_SR_dRcut_highMET_dRl'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_dRl", "hMuTau_SR_dRcut_highMET_dRl;;N", 100, 0, 5)
        h['hMuTau_SR_dRcut_highMET_dRj'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_dRj", "hMuTau_SR_dRcut_highMET_dRj;;N", 100, 0, 5)

        h['hMuTau_SR_dRcut_highMET_lowMt_dRl'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_lowMt_dRl", "hMuTau_SR_dRcut_highMET_lowMt_dRl;;N", 100, 0, 5)
        h['hMuTau_SR_dRcut_highMET_lowMt_dRj'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_lowMt_dRj", "hMuTau_SR_dRcut_highMET_lowMt_dRj;;N", 100, 0, 5)

        h['hMuTau_SR_dRcut_highMET_lowMt_TauPt0'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_lowMt_TauPt0", "hMuTau_SR_dRcut_highMET_lowMt_TauPt0;;N", 500, 0, 500)
        h['hMuTau_SR_dRcut_highMET_lowMt_TauPt1'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_lowMt_TauPt1", "hMuTau_SR_dRcut_highMET_lowMt_TauPt1;;N", 500, 0, 500)
        h['hMuTau_SR_dRcut_highMET_lowMt_TauPt10'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_lowMt_TauPt10", "hMuTau_SR_dRcut_highMET_lowMt_TauPt10;;N", 500, 0, 500)

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
h['hMuCleanedPt_altered'] = ROOT.TH1F ("hMuCleanedPt_altered", "Mu-cleaned Tau P_{T}; P_{T}; a.u.", 500, 0, 500)
h['hECleanedPt'] = ROOT.TH1F ("hECleanedPt", "E-cleaned Tau P_{T}; P_{T}; a.u.", 500, 0, 500)
h['hTauBoostedPt'] = ROOT.TH1F ("hTauBoostedPt", "Boosted Tau P_{T}; P_{T}; a.u.", 500, 0, 500)

h['hMuTau_SR_highMET_lowMt_dRcut_TauPtMass'] = ROOT.TH2F("hMuTau_SR_highMET_lowMt_dRcut_TauPtMass","hMuTau_SR_highMET_lowMt_dRcut_TauPtMass; TauPt(GeV) ; Vis. Mass (GeV)", 500, 0, 500, 100, 0, 100)
h['hMuTau_SR_highMET_lowMt_dRcut_TauPt0Mass'] = ROOT.TH2F("hMuTau_SR_highMET_lowMt_dRcut_TauPt0Mass","hMuTau_SR_highMET_lowMt_dRcut_TauPt0Mass; TauPt(GeV) - DM0 ; Vis. Mass (GeV)", 500, 0, 500, 100, 0, 100)
h['hMuTau_SR_highMET_lowMt_dRcut_TauPt1Mass'] = ROOT.TH2F("hMuTau_SR_highMET_lowMt_dRcut_TauPt1Mass","hMuTau_SR_highMET_lowMt_dRcut_TauPt1Mass; TauPt(GeV) - DM1; Vis. Mass (GeV)", 500, 0, 500, 100, 0, 100)
h['hMuTau_SR_highMET_lowMt_dRcut_TauPt10Mass'] = ROOT.TH2F("hMuTau_SR_highMET_lowMt_dRcut_TauPt10Mass","hMuTau_SR_highMET_lowMt_dRcut_TauPt10Mass; TauPt(GeV) - DM10; Vis. Mass (GeV)", 500, 0, 500, 100, 0, 100)
h['hMuTau_SR_highMET_lowMt_dRcut_NJetMass'] = ROOT.TH2F("hMuTau_SR_highMET_lowMt_dRcut_NJetMass","hMuTau_SR_highMET_lowMt_dRcut_NJetMass; N-Jets ; Vis. Mass (GeV)", 10, 0, 10, 100, 0, 100)

h['hMuTau_SR_highMET_lowMt_dRcut_TauPtMass_loosedR'] = ROOT.TH2F("hMuTau_SR_highMET_lowMt_dRcut_TauPtMass_loosedR","hMuTau_SR_highMET_lowMt_dRcut_TauPtMass_loosedR; TauPt(GeV) ; Vis. Mass_loosedR (GeV)", 500, 0, 500, 100, 0, 100)
h['hMuTau_SR_highMET_lowMt_dRcut_TauPt0Mass_loosedR'] = ROOT.TH2F("hMuTau_SR_highMET_lowMt_dRcut_TauPt0Mass_loosedR","hMuTau_SR_highMET_lowMt_dRcut_TauPt0Mass_loosedR; TauPt(GeV) - DM0 ; Vis. Mass_loosedR (GeV)", 500, 0, 500, 100, 0, 100)
h['hMuTau_SR_highMET_lowMt_dRcut_TauPt1Mass_loosedR'] = ROOT.TH2F("hMuTau_SR_highMET_lowMt_dRcut_TauPt1Mass_loosedR","hMuTau_SR_highMET_lowMt_dRcut_TauPt1Mass_loosedR; TauPt(GeV) - DM1; Vis. Mass_loosedR (GeV)", 500, 0, 500, 100, 0, 100)
h['hMuTau_SR_highMET_lowMt_dRcut_TauPt10Mass_loosedR'] = ROOT.TH2F("hMuTau_SR_highMET_lowMt_dRcut_TauPt10Mass_loosedR","hMuTau_SR_highMET_lowMt_dRcut_TauPt10Mass_loosedR; TauPt(GeV) - DM10; Vis. Mass_loosedR (GeV)", 500, 0, 500, 100, 0, 100)
h['hMuTau_SR_highMET_lowMt_dRcut_NJetMass_loosedR'] = ROOT.TH2F("hMuTau_SR_highMET_lowMt_dRcut_NJetMass_loosedR","hMuTau_SR_highMET_lowMt_dRcut_NJetMass_loosedR; N-Jets ; Vis. Mass_loosedR (GeV)", 10, 0, 10, 100, 0, 100)

h['hMuTau_OS_MetPt_Mt'] = ROOT.TH2F ("hMuTau_OS_MetPt_Mt", "hMuTau_OS_MetPt_Mt ; MET (GeV) ; M_{T} (GeV)", 500, 0, 500, 150, 0, 150)
h['hMuTau_OS_dRcut_MetPt_Mt'] = ROOT.TH2F ("hMuTau_OS_dRcut_MetPt_Mt", "hMuTau_OS_dRcut_MetPt_Mt ; MET (GeV) ; M_{T} (GeV)", 500, 0, 500, 150, 0, 150)
h['hMuTau_SS_MetPt_Mt'] = ROOT.TH2F ("hMuTau_SS_MetPt_Mt", "hMuTau_SS_MetPt_Mt ; MET (GeV) ; M_{T} (GeV)", 500, 0, 500, 150, 0, 150)
h['hMuTau_SS_dRcut_MetPt_Mt'] =ROOT.TH2F ("hMuTau_SS_dRcut_MetPt_Mt", "hMuTau_SS_dRcut_MetPt_Mt ; MET (GeV) ; M_{T} (GeV)", 500, 0, 500, 150, 0, 150)

h['hTauTau_SR_dRcut_highMET_dPhi'] = ROOT.TH2F ("hTauTau_SR_dRcut_highMET_dPhi", "hTauTau_SR_dRcut_highMET_dPhi; dPhi(l,m) ; dPhi(jet,m)", 100, -pi, pi, 100, -pi, pi)
h['hLowTauTau_SR_dRcut_highMET_dPhi'] = ROOT.TH2F ("hLowTauTau_SR_dRcut_highMET_dPhi", "hLowTauTau_SR_dRcut_highMET_dPhi; dPhi(l,m) ; dPhi(jet,m)", 100, -pi, pi, 100, -pi, pi)

book1DHist("hMuTau_SR_dRcut_highMET_lowMt")
book1DHist("hMuTau_SR_dRcut_highMET")
book1DHist("hMuTau_SR_dRcut")
book1DHist("hMuTau_SR")

book1DHist("hMuMu_SR_dRcut_highMET_dPhicut")
book1DHist("hMuMu_SR_dRcut_highMET")
book1DHist("hMuMu_SR_dRcut")

book1DHist("hEMu_SR_dRcut_highMET_dPhicut")

book1DHist("hEMu_SR_dRcut_highMET_dPhicut_allTrig")
book1DHist("hEMu_SR_dRcut_highMET_allTrig")
book1DHist("hEMu_SR_dRcut_allTrig")

book1DHist("hETau_SR_dRcut_highMET_dPhicut")
book1DHist("hETau_SR_dRcut_highMET")
book1DHist("hETau_SR_dRcut")

book1DHist("hLowETau_SR_dRcut_highMET_dPhicut")
book1DHist("hLowETau_SR_dRcut_highMET")
book1DHist("hLowETau_SR_dRcut")

book1DHist("hTauTau_SR_dRcut_highMET")
book1DHist("hTauTau_SR_dRcut")

book1DHist("hLowTauTau_SR_dRcut_highMET")
book1DHist("hLowTauTau_SR_dRcut")

book1DHist("hTauTau_dRcut_lowMET")
book1DHist("hTauTau_dRcut")

book1DHist("hLowTauTau_dRcut_lowMET")
book1DHist("hLowTauTau_dRcut")

book1DHist("hEE_SR_dRcut")
book1DHist("hEE_SR_dRcut_highMET")
book1DHist("hEE_SR_dRcut_highMET_dPhicut")


for key in h.keys():
    h[key].Sumw2()


def EE_Channel(ele):

    isEE = 0 

    e1 = ROOT.TLorentzVector()
    e1.SetPtEtaPhiM(ele[0].pt, ele[0].eta, ele[0].phi, ele[0].mass)
    e2 = ROOT.TLorentzVector()
    e2.SetPtEtaPhiM(ele[1].pt, ele[1].eta, ele[1].phi, ele[1].mass)

    j = ROOT.TLorentzVector()
    j.SetPtEtaPhiM(s_j[0].pt, s_j[0].eta, s_j[0].phi, s_j[0].mass)

    m = ROOT.TLorentzVector()
    m.SetPtEtaPhiM(met_pt, 0, met_phi, 0)

    if (e1+e2).M() > 1.0:
        if e1.DeltaR(e2) < 0.4 and e1.DeltaR(j)> 0.8  and e2.DeltaR(j) > 0.8:
            plot1Dhist('hEE_SR_dRcut', e1, e2, j, m)
            if met_pt > 100:
                plot1Dhist('hEE_SR_dRcut_highMET', e1, e2, j, m)
                if abs(m.DeltaPhi(e1)) < 1.0 and abs(m.DeltaPhi(j)) > 2.0:
                    if ( ( j.Pt() > 500 and ( isHT == 1 ) ) or ( e1.Pt() > 35 and isIsoEle == 1 ) ) and isData == 0 :
                        isEE = 1
                        plot1Dhist('hEE_SR_dRcut_highMET_dPhicut', e1, e2, j, m)
                        h['hTauEvents'].Fill(1,1)
                        h['ChOverlap'].Fill(9,1)
                        if isEMu == 1 or isMuMu == 1: h['ChOverlap'].Fill(3,1)


    return isEE

def MuMu_Channel(mu):

   isMuMu = 0
   h['hMuMu_Events'].Fill(1, genweight)

   mu1 = ROOT.TLorentzVector()
   mu1.SetPtEtaPhiM(mu[0].pt, mu[0].eta, mu[0].phi, mu[0].mass)

   mu2 = ROOT.TLorentzVector()
   mu2.SetPtEtaPhiM(mu[1].pt, mu[1].eta, mu[1].phi, mu[1].mass)

   j = ROOT.TLorentzVector()
   j.SetPtEtaPhiM(s_j[0].pt, s_j[0].eta, s_j[0].phi, s_j[0].mass)

   m = ROOT.TLorentzVector()
   m.SetPtEtaPhiM(met_pt, 0, met_phi, 0)

   if j.Pt() > 100 and (mu1+mu2).M() > 1:   

       isJetHTEvent = 0
       isSingleMuonEvent = 0

       if met_pt < 100 and len(s_b) == 0:
           if j.Pt() > 500 and isSingleJet == 1 :
               h['hMuMu_Trigger'].Fill(1, genweight)
           if j.Pt() > 500 and isHT == 1 :
               h['hMuMu_Trigger'].Fill(2, genweight)
           if j.Pt() > 600 and isSingleJet == 1 :
               h['hMuMu_Trigger'].Fill(3, genweight)
           if j.Pt() > 600 and isHT == 1 :
               h['hMuMu_Trigger'].Fill(4, genweight)
           if ( j.Pt() > 500 and ( isSingleJet == 1 or isHT == 1 ) ) :
               h['hMuMu_Trigger'].Fill(5, genweight)
           if ( j.Pt() > 600 and ( isSingleJet == 1 or isHT == 1 ) ) :
               h['hMuMu_Trigger'].Fill(6, genweight)

       if plotHTTrig == 1 :
           if ( j.Pt() > 500 and isHT == 1 ) : isJetHTEvent = 1

       if plotJetHT == 1 :
           if ( j.Pt() > 600 and ( isSingleJet == 1 or isHT == 1 ) ) : isJetHTEvent = 1

       if plotSingleMuon == 1 :
           if ( mu1.Pt() > 50 and isMu == 1 ) or ( mu1.Pt() > 27 and isIsoMu == 1 ) : isSingleMuonEvent = 1

       if met_pt < 100 and len(s_b) == 0:
           if isJetHTEvent == 1 and isSingleMuonEvent == 1:
               h['hMuMu_Inclusive_Trigger'].Fill(1, genweight)

           if isJetHTEvent == 1 and isSingleMuonEvent == 0:
               h['hMuMu_Inclusive_Trigger'].Fill(2, genweight)

           if isJetHTEvent == 0 and isSingleMuonEvent == 1:
               h['hMuMu_Inclusive_Trigger'].Fill(3, genweight)

           if isJetHTEvent == 1 or isSingleMuonEvent == 1:
               h['hMuMu_Inclusive_Trigger'].Fill(4, genweight)

       if ( isData == 0 and ( ( plotHTTrig == 1 and isJetHTEvent == 1 ) or ( plotJetHT == 1 and isJetHTEvent == 1 ) or ( plotSingleMuon == 1 and isSingleMuonEvent == 1 ) ) ) or \
          ( isData == 1 and ( isJetHTDataset == 1 and ( isJetHTEvent == 1 and isSingleMuonEvent == 0 ) ) ) or \
          ( isData == 1 and ( isSingleMuonDataset == 1 and isSingleMuonEvent == 1 ) ):

           # print("data? :", isData)
           # print("isJetHTDataset, isJetHTEvent, isSingleMuonEvent :", isJetHTDataset, isJetHTEvent, isSingleMuonEvent)
           # print("isSingleMuonDataset, isJetHTEvent, isSingleMuonEvent, ", isSingleMuonDataset, isJetHTEvent , isSingleMuonEvent)

           h['hMuMu_Events'].Fill(2, genweight)

           if len(s_b) == 0:
                if mu1.DeltaR(mu2) < 0.4 and mu1.DeltaR(j) > 0.8  and mu2.DeltaR(j) > 0.8:
                   if isData == 0 :
                       plot1Dhist('hMuMu_SR_dRcut', mu1, mu2, j, m)
                   if met_pt > 100:
                       if isData == 0 :
                           plot1Dhist('hMuMu_SR_dRcut_highMET', mu1, mu2, j, m)
                       if abs(m.DeltaPhi(mu1)) < 1 and abs(m.DeltaPhi(j)) > 2:
                           isMuMu = 1
                           h['ChOverlap'].Fill(7,1)
                           if isData == 0 :
                               h['hTauEvents'].Fill(2,genweight)
                               plot1Dhist('hMuMu_SR_dRcut_highMET_dPhicut', mu1, mu2, j, m)
#                           h['hMuMu_SR_dRcut_highMET_dPhicut'].Fill((mu1+mu2).M(), genweight)


   return isMuMu

def EMu_Channel(ele,mu_emu):

   isEMu = 0
   h['hEMu_Events'].Fill(1, genweight)

   mu = ROOT.TLorentzVector()
   mu.SetPtEtaPhiM(mu_emu[0].pt, mu_emu[0].eta, mu_emu[0].phi, mu_emu[0].mass)

   e = ROOT.TLorentzVector()
   e.SetPtEtaPhiM(ele[0].pt, ele[0].eta, ele[0].phi, ele[0].mass)

   j = ROOT.TLorentzVector()
   j.SetPtEtaPhiM(s_j[0].pt, s_j[0].eta, s_j[0].phi, s_j[0].mass)

   m = ROOT.TLorentzVector()
   m.SetPtEtaPhiM(met_pt, 0, met_phi, 0)

   isJetHTEvent = 0
   isSingleMuonEvent = 0

   if met_pt < 100 and (e+mu).M() > 1:
       if j.Pt() > 500 and isSingleJet == 1 :
           h['hEMu_Trigger'].Fill(1, genweight)
       if j.Pt() > 500 and isHT == 1 :
           h['hEMu_Trigger'].Fill(2, genweight)
       if j.Pt() > 600 and isSingleJet == 1 :
           h['hEMu_Trigger'].Fill(3, genweight)
       if j.Pt() > 600 and isHT == 1 :
           h['hEMu_Trigger'].Fill(4, genweight)
       if ( j.Pt() > 500 and ( isSingleJet == 1 or isHT == 1 ) ) :
           h['hEMu_Trigger'].Fill(5, genweight)
       if ( j.Pt() > 600 and ( isSingleJet == 1 or isHT == 1 ) ) :
           h['hEMu_Trigger'].Fill(6, genweight)

   if plotHTTrig == 1 and (e+mu).M() > 1 :
       if ( j.Pt() > 500 and isHT == 1 ) : isJetHTEvent = 1

   if plotJetHT == 1 and (e+mu).M() > 1:
       if ( j.Pt() > 600 and ( isSingleJet == 1 or isHT == 1 ) ) : isJetHTEvent = 1

   if ( mu.Pt() > 50 and isMu == 1 ) or ( mu.Pt() > 27 and isIsoMu == 1 ) : isSingleMuonEvent = 1

   if met_pt < 100 and len(s_b) == 0 and j.Pt() > 100 and (e+mu).M() > 1:
       if isJetHTEvent == 1 and isSingleMuonEvent == 1:
           h['hEMu_Inclusive_Trigger'].Fill(1, genweight)

       if isJetHTEvent == 1 and isSingleMuonEvent == 0:
           h['hEMu_Inclusive_Trigger'].Fill(2, genweight)

       if isJetHTEvent == 0 and isSingleMuonEvent == 1:
           h['hEMu_Inclusive_Trigger'].Fill(3, genweight)

       if isJetHTEvent == 1 or isSingleMuonEvent == 1:
           h['hEMu_Inclusive_Trigger'].Fill(4, genweight)


   if isData == 0 and j.Pt() > 100 :

       if ( j.Pt() > 500 and ( isHT == 1 ) ) \
          or ( isMuonEG == 1 and ( ( mu.Pt() > 8 and e.Pt() > 23 ) or ( mu.Pt() > 23 and e.Pt() > 12 ) ) ) \
          or ( mu.Pt() > 50 and isMu == 1 ) or ( mu.Pt() > 27 and isIsoMu == 1 ) \
          or ( e.Pt() > 35 and (isIsoEle == 1) ) :

           if (e+mu).M() > 1.0 :
               if mu.DeltaR(e) < 0.4 and mu.DeltaR(j)> 0.8  and e.DeltaR(j) > 0.8 :
                   plot1Dhist('hEMu_SR_dRcut_allTrig',e,mu,j,m)
                   if met_pt > 100.0 :
                       plot1Dhist('hEMu_SR_dRcut_highMET_allTrig',e,mu,j,m)
                       if abs(m.DeltaPhi(mu)) < 1.0 and abs( m.DeltaPhi(j) ) > 2.0 :
                           plot1Dhist('hEMu_SR_dRcut_highMET_dPhicut_allTrig',e,mu,j,m)
                           isEMu = 1

                    
   return isEMu


def Mt(lepton, met):

    cos = np.cos(met.DeltaPhi(lepton))
    Mt = np.sqrt(2*lepton.Pt()*met.Pt()*(1-cos))

    return Mt

def MuTau_Channel(mtau, lmuon):

       isMuTau = 0
       h['hMuTau_Events'].Fill(1, genweight)

       mu = ROOT.TLorentzVector()
       mu.SetPtEtaPhiM(lmuon[0].pt, lmuon[0].eta, lmuon[0].phi, lmuon[0].mass)

       tau = ROOT.TLorentzVector()
       tau.SetPtEtaPhiM(mtau[0].pt, mtau[0].eta, mtau[0].phi, mtau[0].mass)

       if len(s_mu) > 1:
           mu2 = ROOT.TLorentzVector()
           mu2.SetPtEtaPhiM(lmuon[1].pt, lmuon[1].eta, lmuon[1].phi, lmuon[1].mass)

       if len(s_j) > 1 :
           j0 = ROOT.TLorentzVector()
           j0.SetPtEtaPhiM(s_j[0].pt, s_j[0].eta, s_j[0].phi, s_j[0].mass)
           if tau.DeltaR(j0) < 0.2 :
               j = ROOT.TLorentzVector()
               j.SetPtEtaPhiM(s_j[1].pt, s_j[1].eta, s_j[1].phi, s_j[1].mass)
               j2 = ROOT.TLorentzVector()
               j2.SetPtEtaPhiM(s_j[0].pt, s_j[0].eta, s_j[0].phi, s_j[0].mass)
           else:
               j = j0
               j2 = ROOT.TLorentzVector()
               j2.SetPtEtaPhiM(s_j[1].pt, s_j[1].eta, s_j[1].phi, s_j[1].mass)
       if len(s_j) == 1:
           j = ROOT.TLorentzVector()
           j.SetPtEtaPhiM(s_j[0].pt, s_j[0].eta, s_j[0].phi, s_j[0].mass)

       m = ROOT.TLorentzVector()
       m.SetPtEtaPhiM(met_pt, 0, met_phi, 0)

       if lmuon[0].charge*mtau[0].charge < 0 and len(s_b) == 0 and met_pt < 100 and (mu+tau).M() > 1:
           if j.Pt() > 500 and isSingleJet == 1 :
               h['hMuTau_Trigger'].Fill(1, genweight)
           if j.Pt() > 500 and isHT == 1 :
               h['hMuTau_Trigger'].Fill(2, genweight)
           if j.Pt() > 600 and isSingleJet == 1 :
               h['hMuTau_Trigger'].Fill(3, genweight)
           if j.Pt() > 600 and isHT == 1 :
               h['hMuTau_Trigger'].Fill(4, genweight)
           if ( j.Pt() > 500 and ( isSingleJet == 1 or isHT == 1 ) ) :
               h['hMuTau_Trigger'].Fill(5, genweight)
           if ( j.Pt() > 600 and ( isSingleJet == 1 or isHT == 1 ) ) :
               h['hMuTau_Trigger'].Fill(6, genweight)

       isJetHTEvent = 0
       isSingleMuonEvent = 0

       if plotHTTrig == 1 and (mu+tau).M() > 1:
           if ( j.Pt() > 500 and isHT == 1 ) : isJetHTEvent = 1

       if plotJetHT == 1 and (mu+tau).M() > 1:
           if ( j.Pt() > 600 and ( isSingleJet == 1 or isHT == 1 ) ) : isJetHTEvent = 1

       if ( mu.Pt() > 50 and isMu == 1 ) or ( mu.Pt() > 27 and isIsoMu == 1 ) : isSingleMuonEvent = 1

       if lmuon[0].charge*mtau[0].charge < 0 and len(s_b) == 0 and met_pt < 100 and j.Pt() > 100 and (mu+tau).M() > 1:
           if isJetHTEvent == 1 and isSingleMuonEvent == 1:
               h['hMuTau_Inclusive_Trigger'].Fill(1, genweight)
           if isJetHTEvent == 1 and isSingleMuonEvent == 0:
               h['hMuTau_Inclusive_Trigger'].Fill(2, genweight)
           if isJetHTEvent == 0 and isSingleMuonEvent == 1:
               h['hMuTau_Inclusive_Trigger'].Fill(3, genweight)
           if isJetHTEvent == 1 or isSingleMuonEvent == 1:
               h['hMuTau_Inclusive_Trigger'].Fill(4, genweight)

       if ( isData == 0 and ( ( plotHTTrig == 1 and isJetHTEvent == 1 ) or ( plotJetHT == 1 and isJetHTEvent == 1 ) or ( plotSingleMuon == 1 and isSingleMuonEvent == 1 ) ) ) or \
          ( isData == 1 and ( isJetHTDataset == 1 and ( isJetHTEvent == 1 and isSingleMuonEvent == 0 ) ) ) or \
          ( isData == 1 and ( isSingleMuonDataset == 1 and isSingleMuonEvent == 1 ) ):

           if isData == 0 or isAltered == 1:
                if j.Pt() > 100 and (mu+tau).M() > 1:
                    if lmuon[0].charge*mtau[0].charge < 0 and len(s_b) == 0 :
                        h['hMuTau_SR_Events'].Fill(1, genweight)
                        h['hMuTau_SR_dRl'].Fill(mu.DeltaR(tau), genweight)
                        h['hMuTau_SR_dRj'].Fill(j.DeltaR(mu+tau), genweight)
                        plot1Dhist("hMuTau_SR", tau, mu, j, m)

                        if mu.DeltaR(tau) < 0.4 and mu.DeltaR(j) > 0.8 and tau.DeltaR(j) > 0.8 and mu.DeltaR(tau) > 0.05 :
                            h['hMuTau_SR_Events'].Fill(2, genweight)
                            h['hMuTau_SR_dRcut_dRl'].Fill(mu.DeltaR(tau), genweight)
                            h['hMuTau_SR_dRcut_dRj'].Fill(j.DeltaR(mu+tau), genweight)
                            plot1Dhist("hMuTau_SR_dRcut", tau, mu, j, m)

                            if met_pt > 100 :
                                h['hMuTau_SR_Events'].Fill(3, genweight)
                                h['hMuTau_SR_dRcut_highMET_dRl'].Fill(mu.DeltaR(tau), genweight)
                                h['hMuTau_SR_dRcut_highMET_dRj'].Fill(j.DeltaR(mu+tau), genweight)
                                plot1Dhist("hMuTau_SR_dRcut_highMET", tau, mu, j, m)

                                if Mt(mu,m) < 50 :
                                    h['hMuTau_SR_Events'].Fill(4, genweight)
                                    h['hMuTau_SR_dRcut_highMET_lowMt_dRl'].Fill(mu.DeltaR(tau), genweight)
                                    h['hMuTau_SR_dRcut_highMET_lowMt_dRj'].Fill(j.DeltaR(mu+tau), genweight)
                                    plot1Dhist("hMuTau_SR_dRcut_highMET_lowMt", tau, mu, j, m)
                                    isMuTau = 1
                                    h['ChOverlap'].Fill(10,1)
                                    if isEMu == 1 or isMuMu == 1 or isEE == 1: h['ChOverlap'].Fill(4,1)
                                    h['hTauEvents'].Fill(4,genweight)

                                    if mtau[0].decaymode == 0 : h['hMuTau_SR_dRcut_highMET_lowMt_TauPt0'].Fill(tau.Pt(), genweight)
                                    if mtau[0].decaymode == 1 : h['hMuTau_SR_dRcut_highMET_lowMt_TauPt1'].Fill(tau.Pt(), genweight)
                                    if mtau[0].decaymode == 10 : h['hMuTau_SR_dRcut_highMET_lowMt_TauPt10'].Fill(tau.Pt(), genweight)

                                    if isJetHTEvent == 1 : h['hMuTau_SR_Trigger'].Fill(7, genweight)
                                    h['hMuTau_SR_Inclusive_Trigger'].Fill(7, genweight)

                                    h["hMuTau_SR_highMET_lowMt_dRcut_TauPtMass"].Fill(tau.Pt(), (mu+tau).M(), genweight)
                                    h["hMuTau_SR_highMET_lowMt_dRcut_NJetMass"].Fill(len(s_j), (mu+tau).M(), genweight)

                                    if mtau[0].decaymode == 0 : h["hMuTau_SR_highMET_lowMt_dRcut_TauPt0Mass"].Fill(tau.Pt(), (mu+tau).M(), genweight)
                                    if mtau[0].decaymode == 1 : h["hMuTau_SR_highMET_lowMt_dRcut_TauPt1Mass"].Fill(tau.Pt(), (mu+tau).M(), genweight)
                                    if mtau[0].decaymode == 10 : h["hMuTau_SR_highMET_lowMt_dRcut_TauPt10Mass"].Fill(tau.Pt(), (mu+tau).M(), genweight)

           if j.Pt() > 100 and (mu+tau).M() > 1:
                if lmuon[0].charge*mtau[0].charge < 0 and len(s_b) == 0 :
                    if mu.DeltaR(tau) < 0.4 and mu.DeltaR(j) > 0.8 and tau.DeltaR(j) > 0.8 and mu.DeltaR(tau) > 0.05:
                        if met_pt > 100 :
                            if Mt(mu,m) < 50 :
                                if ( j.Pt() > 500 and isSingleJet == 1 ) or isSingleMuonEvent == 1 :
                                    h['hMuTau_SR_Inclusive_Trigger'].Fill(1, genweight)
                                if ( j.Pt() > 500 and isHT == 1 ) or isSingleMuonEvent == 1 :
                                    h['hMuTau_SR_Inclusive_Trigger'].Fill(2, genweight)
                                if ( j.Pt() > 600 and isSingleJet == 1 ) or isSingleMuonEvent == 1 :
                                    h['hMuTau_SR_Inclusive_Trigger'].Fill(3, genweight)
                                if ( j.Pt() > 600 and isHT == 1 ) or isSingleMuonEvent == 1 :
                                    h['hMuTau_SR_Inclusive_Trigger'].Fill(4, genweight)
                                if ( ( j.Pt() > 500 and ( isSingleJet == 1 or isHT == 1 ) ) ) or isSingleMuonEvent == 1 :
                                    h['hMuTau_SR_Inclusive_Trigger'].Fill(5, genweight)
                                if ( ( j.Pt() > 600 and ( isSingleJet == 1 or isHT == 1 ) ) ) or isSingleMuonEvent == 1 :
                                    h['hMuTau_SR_Inclusive_Trigger'].Fill(6, genweight)

                                if ( j.Pt() > 500 and isSingleJet == 1 ) : 
                                    h['hMuTau_SR_Trigger'].Fill(1, genweight)
                                if ( j.Pt() > 500 and isHT == 1 ) :
                                    h['hMuTau_SR_Trigger'].Fill(2, genweight)
                                if ( j.Pt() > 600 and isSingleJet == 1 ) :
                                    h['hMuTau_SR_Trigger'].Fill(3, genweight)
                                if ( j.Pt() > 600 and isHT == 1 ) :
                                    h['hMuTau_SR_Trigger'].Fill(4, genweight)
                                if ( ( j.Pt() > 500 and ( isSingleJet == 1 or isHT == 1 ) ) ) :
                                    h['hMuTau_SR_Trigger'].Fill(5, genweight)
                                if ( ( j.Pt() > 600 and ( isSingleJet == 1 or isHT == 1 ) ) ) :
                                    h['hMuTau_SR_Trigger'].Fill(6, genweight)


           # if plot2DforTau == 1:
           #     if j.Pt() > 1 and (mu+tau).M() > 1:
           #         if lmuon[0].charge*mtau[0].charge < 0 and len(s_b) == 0 :
           #             plot2DTauPt("hMuTau_SS")


       return isMuTau


def plot1Dhist(region,l1,l2,j,m):

    h[region+"_Lepton1Pt"].Fill(l1.Pt(), genweight)
    h[region+"_Mass"].Fill((l1+l2).M(), genweight)
    h[region+"_Lepton2Pt"].Fill(l2.Pt(), genweight)
    h[region+"_JetPt"].Fill(j.Pt(), genweight)
    h[region+"_Mt"].Fill(Mt(l2,m), genweight)
    h[region+"_MetPt"].Fill(m.Pt(), genweight)
    h[region+"_Nj"].Fill(len(s_j), genweight)


def ETau_Channel(etau, lelectron, reco):

   isETauTight = 0

   e = ROOT.TLorentzVector()
   e.SetPtEtaPhiM(lelectron[0].pt, lelectron[0].eta, lelectron[0].phi, lelectron[0].mass)

   tau = ROOT.TLorentzVector()
   tau.SetPtEtaPhiM(etau[0].pt, etau[0].eta, etau[0].phi, etau[0].mass)

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

   if len(s_b) == 0 and etau[0].charge*lelectron[0].charge < 0 and isData == 0 and j.Pt() > 100 :

       if ( j.Pt() > 500 and ( isHT == 1 ) ) \
          or ( e.Pt() > 35 and (isIsoEle == 1 ) ) \
          or ( j.Pt() > 500 and met_pt > 200 and isHTMHT == 1 ) :
           
           if e.DeltaR(tau) < 0.4 and e.DeltaR(j) > 0.8 and tau.DeltaR(j) > 0.8 :
               if e.DeltaR(tau) > 0.05 :
                   if reco == "standard" : plot1Dhist("hETau_SR_dRcut", tau, e, j, m)
                   if reco == "lowPt" : plot1Dhist("hLowETau_SR_dRcut", tau, e, j, m)
                   if met_pt > 100:
                       if reco == "standard" : plot1Dhist("hETau_SR_dRcut_highMET", tau, e, j, m)
                       if reco == "lowPt" : plot1Dhist("hLowETau_SR_dRcut_highMET", tau, e, j, m)
                       if abs(m.DeltaPhi(e)) < 1 and abs(m.DeltaPhi(j)) > 2:
                           isETauTight = 1
                           h['ChOverlap'].Fill(11,1)
                           if isMuTau == 1 or isEMu == 1 or isMuMu == 1 or isEE == 1: h['ChOverlap'].Fill(5,1)
                           if reco == "standard" : 
                               plot1Dhist("hETau_SR_dRcut_highMET_dPhicut", tau, e, j, m)
                               h['hTauEvents'].Fill(5,genweight)
                           if reco == "lowPt" : 
                               plot1Dhist("hLowETau_SR_dRcut_highMET_dPhicut", tau, e, j, m)
                               h['hTauEvents'].Fill(6,genweight)

   return isETauTight


def TauTau_Channel(boosted, reco):

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

   isJetHTEvent = 0

   if plotHTTrig == 1 and (tau1+tau2).M() > 1:
       if ( j.Pt() > 500 and isHT == 1 ) : isJetHTEvent = 1

   if (tau1+tau2).M() > 1 and isData == 0 and j.Pt() > 100 :
       if ( j.Pt() > 500 and isHT == 1 ) \
          or ( j.Pt() > 500 and met_pt > 200 and isHTMHT == 1 ) :

           if tau1.DeltaR(tau2) < 0.4 and tau1.DeltaR(j) > 0.8 and tau2.DeltaR(j) > 0.8:
               if reco == "standard" : plot1Dhist("hTauTau_SR_dRcut", tau1, tau2, j, m)
               if reco == "lowPt" : plot1Dhist("hLowTauTau_SR_dRcut", tau1, tau2, j, m)

               if met_pt > 100:
                   if reco == "standard" : plot1Dhist("hTauTau_SR_dRcut_highMET", tau1, tau2, j, m)
                   if reco == "lowPt" : plot1Dhist("hLowTauTau_SR_dRcut_highMET", tau1, tau2, j, m)

                   h['ChOverlap'].Fill(12,1)
                   if isETauTight == 1 or isMuMu == 1 or isEMu == 1 or isMuTau == 1 or isEE == 1 : h['ChOverlap'].Fill(6,1)
                   if reco == "standard" :
                       h['hTauEvents'].Fill(7,genweight)
                       h["hTauTau_SR_dRcut_highMET_dPhi"].Fill(tau1.DeltaPhi(m),j.DeltaPhi(m), genweight)
                   if reco == "standard" :
                       h['hTauEvents'].Fill(9,genweight)
                       h["hLowTauTau_SR_dRcut_highMET_dPhi"].Fill(tau1.DeltaPhi(m),j.DeltaPhi(m), genweight)



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
lowPtElectrons = ROOT.ElectronInfoDS()
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
fchain.SetBranchAddress("LowPtElectrons", ROOT.AddressOf(lowPtElectrons))

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

   genweight = 1

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

   unclean = []
   eclean = []
   mclean = []
   lowEclean = []

   boosted = []
   mclean_altered = []

   if tausUnCleaned.size()>0:
       for i in range(tausUnCleaned.size()):
           tau = tausUnCleaned.at(i)
           if abs(tau.eta) < 2.3 :
               if tau.mvaid >= 4:
                   h['hTauUnCleanedPt'].Fill(tau.pt, genweight)
                   unclean+=[tau]

   if tausECleaned.size()>0:
       for i in range(tausECleaned.size()):
           tau = tausECleaned.at(i)
           if abs(tau.eta) < 2.3 :
               if tau.mvaid >= 4:
                   h['hECleanedPt'].Fill(tau.pt, genweight)
                   eclean+=[tau]

   if tausMCleaned.size()>0:
       for i in range(tausMCleaned.size()):
           tau = tausMCleaned.at(i)
           if abs(tau.eta) < 2.3 :
               if tau.mvaid >= 1 and tau.mvaid < 4 :
                   mclean_altered+=[tau]
                   h['hMuCleanedPt_altered'].Fill(tau.pt, genweight)
               if tau.mvaid >= 4 :
                   h['hMuCleanedPt'].Fill(tau.pt, genweight)
                   mclean+=[tau]

   if tausLowPtECleaned.size()>0:
       for i in range(tausLowPtECleaned.size()):
           tau = tausLowPtECleaned.at(i)
           if tau.mvaid >= 4:
               lowEclean+=[tau]

   if tausBoosted.size()>0:
       for i in range(tausBoosted.size()):
           tau = tausBoosted.at(i)
           if abs(tau.eta) < 2.3 :
               if tau.mvaid >= 1:
                   h['hTauBoostedPt'].Fill(tau.pt, genweight)
                   boosted+=[tau]

   if jets.size()>0:
      for i in range(jets.size()):
         jet = jets.at(i)
         if jet.id == 2:
             h['hJetPt'].Fill(jet.pt, genweight)
             h['hBtag'].Fill(jet.deepjet)
             s_j+=[jet]
             if jet.deepjet > 0.7476:
                 h['hBJetPt'].Fill(jet.pt, genweight)
                 s_b+=[jet]

   if muons.size()>0:
      for i in range(muons.size()):
         muon = muons.at(i)
         if abs(muon.eta) < 2.4 : 
             if muon.id >= 1:
                 h['hMuonPt'].Fill(muon.pt, genweight) 
                 s_mu+=[muon]
                 if muon.iso < 0.25:
                     h['hIsoMuonPt'].Fill(muon.pt, genweight)
                     s_isomu+=[muon]

   if electrons.size()>0:
      for i in range(electrons.size()):
         electron = electrons.at(i)
         if abs(electron.eta) < 2.5 :
             if electron.id >= 1 :
                 h['hElectronPt'].Fill(electron.pt, genweight)
                 s_e+=[electron]
                 if electron.iso >= 1:
                     h['hIsoElectronPt'].Fill(electron.pt, genweight)
                     s_isoe+=[electron]

   if lowPtElectrons.size() > 0 :
       for i in range(lowPtElectrons.size()):
           lowPtElectron = lowPtElectrons.at(i)
           if abs(lowPtElectron.eta) < 2.5 :
               s_lowE += [lowPtElectron]

   s_j.sort(key=lambda x: x.pt, reverse=True)
   s_e.sort(key=lambda x: x.pt, reverse=True)
   s_lowE.sort(key=lambda x: x.pt, reverse=True)
   s_isoe.sort(key=lambda x: x.pt, reverse=True)
   s_mu.sort(key=lambda x: x.pt, reverse=True)
   s_isomu.sort(key=lambda x: x.pt, reverse=True)
   unclean.sort(key=lambda x: x.pt, reverse=True)
   eclean.sort(key=lambda x: x.pt, reverse=True)
   mclean.sort(key=lambda x: x.pt, reverse=True)
   mclean_altered.sort(key=lambda x: x.pt, reverse=True)
   lowEclean.sort(key=lambda x: x.pt, reverse=True)

   boosted.sort(key=lambda x: x.pt, reverse=True)

   isEMu = 0
   isMuMu = 0
   isEE = 0
   isMuTau = 0
   isETauTight = 0
   isLowETauTight = 0

   if len(s_isomu) > 1 and len(s_j) > 0 and s_isomu[0].charge*s_isomu[1].charge < 0 : 
       if MuMu_Channel(s_isomu) == 1: continue

   if len(s_isomu) > 0 and len(s_isoe) > 0 and len(s_j) > 0 and len(s_b) == 0 and s_isoe[0].charge*s_isomu[0].charge < 0 : 
       if EMu_Channel(s_isoe,s_isomu) == 1: continue

   if len(s_mu) > 0 and len(mclean) > 0 and len(s_j) > 0 : 
       if MuTau_Channel(mclean, s_mu) == 1 : continue

   if len(s_isoe) > 1 and len(s_j) > 0 and len(s_b) == 0 and s_isoe[0].charge*s_isoe[1].charge < 0 :                                                     
       if EE_Channel(s_isoe) == 1: continue      

   if len(s_e) > 0 and len(eclean) > 0 and len(s_j) > 0 and len(s_b) == 0 and eclean[0].charge*s_e[0].charge < 0:
       if ETau_Channel(eclean, s_e, "standard") == 0 : 
           if len(boosted) > 1 and len(s_j) > 0 and len(s_b) == 0 and boosted[0].charge*boosted[1].charge < 0 :
               TauTau_Channel(boosted, "standard")

   if len(s_lowE) > 0 and len(lowEclean) > 0 and len(s_j) > 0 and len(s_b) == 0 and lowEclean[0].charge*s_lowE[0].charge < 0:
       if ETau_Channel(lowEclean, s_lowE, "lowPt") == 0 :
           if len(boosted) > 1 and len(s_j) > 0 and len(s_b) == 0 and boosted[0].charge*boosted[1].charge < 0 :
               TauTau_Channel(boosted, "lowPt")

           

   # if isAltered == 1 and plot2DforTau == 1:
   #     if len(s_mu) > 0 and len(mclean_altered) > 0 and len(s_j) > 0:
   #         plot2DTauPt(mclean_altered, s_mu)

   # if isAltered == 0 and plot2DforTau == 1:
   #     if len(s_mu) > 0 and len(mclean) > 0 and len(s_j) > 0:
   #         plot2DTauPt(mclean, s_mu)


#   if len(s_isoe) > 1 and len(s_j) > 0 and len(s_b) == 0 and s_isoe[0].charge*s_isoe[1].charge < 0 : 
#       if EE_Channel(s_isoe) == 1: continue



out.cd()

for key in h.keys():
    h[key].Write()

out.Close()

print("--- %s seconds ---" % (time.time() - start_time))
