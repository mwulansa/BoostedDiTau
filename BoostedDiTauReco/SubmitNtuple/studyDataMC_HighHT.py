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

plot2DforTau = 1

isAltered = 1
print("isAltered = ", isAltered)

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

if plotSingleMuon == 1:
    outputFileName = outputFileDir+"h_debugMuTau_HighHT_"+inputFileListName.split("/")[-1].replace(".txt",".root")

if plotJetHT == 1:
    outputFileName = outputFileDir+"h_debugMuTau_HighHT_JetHT_"+inputFileListName.split("/")[-1].replace(".txt",".root")

if isData == 1 :
    outputFileName = outputFileDir+"h_debugMuTau_HighHT_Data_"+inputFileListName.split("/")[-1].replace(".txt",".root")
    if isAltered == 1 :
        outputFileName = outputFileDir+"h_debugMuTau_HighHT_Data_Altered_"+inputFileListName.split("/")[-1].replace(".txt",".root")
    if plotSemiLeptonic == 1:
        outputFileName = outputFileDir+"h_debugMuTau_HighHT_FullyLeptonic_Data_"+inputFileListName.split("/")[-1].replace(".txt",".root")
        if isAltered == 1 :
            outputFileName = outputFileDir+"h_debugMuTau_HighHT_FullyLeptonic_Data_Altered_"+inputFileListName.split("/")[-1].replace(".txt",".root")

if plotSemiLeptonic == 1 and plotJetHT == 1 and isData == 0 :
    outputFileName = outputFileDir+"h_debugMuTau_HighHT_FullyLeptonic_JetHT_"+inputFileListName.split("/")[-1].replace(".txt",".root")

if plotSemiLeptonic == 1 and plotHTTrig == 1 and isData == 0 :
    outputFileName = outputFileDir+"h_debugMuTau_HighHT_FullyLeptonic_HTTrig_"+inputFileListName.split("/")[-1].replace(".txt",".root")

if plotSingleMuonHTTrig == 1 and isData== 0 :
    if plotSemiLeptonic == 1:
        outputFileName = outputFileDir+"h_debugMuTau_HighHT_FullyLeptonic_Inclusive_"+inputFileListName.split("/")[-1].replace(".txt",".root")
        if isAltered == 1 :
            outputFileName = outputFileDir+"h_debugMuTau_HighHT_FullyLeptonic_Inclusive_Altered_"+inputFileListName.split("/")[-1].replace(".txt",".root")
    if plotSemiLeptonic == 0:
        outputFileName = outputFileDir+"h_debugMuTau_HighHT_Inclusive_"+inputFileListName.split("/")[-1].replace(".txt",".root")
        if isAltered == 1 :
            outputFileName = outputFileDir+"h_debugMuTau_HighHT_Inclusive_Altered_"+inputFileListName.split("/")[-1].replace(".txt",".root")


out=ROOT.TFile.Open(outputFileName,'recreate')
print(outputFileName)

fchain = ROOT.TChain('tcpNtuples/analysisTree')
chain2 = ROOT.TChain('tcpTrigNtuples/triggerTree')
if isData == 0:
    chain3 = ROOT.TChain('lumiSummary/lumiTree')
#gchain = ROOT.TChain('tcpGenNtuples/genTree')

pi = np.pi

h = {}

h['hEvents'] = ROOT.TH1F ("NEvents", "Number of Events; ;N", 2, 0, 2)
h['hTauEvents'] = ROOT.TH1F ("hTauEvents", "Number of each decay modes;;a.u.",7,0,7)

h['hMuMu_Events'] = ROOT.TH1F ("hMuMu_Events", "hMuMu_Events;;N", 6, 1, 7)
h['hEMu_Events'] = ROOT.TH1F ("hEMu_Events", "hEMu_Events;;N", 6, 1, 7)
h['hMuTau_Events'] = ROOT.TH1F ("hMuTau_Events", "hMuTau_Events;;N", 6, 1, 7)
h['hMuTau_SR_Events'] = ROOT.TH1F ("hMuTau_SR_Events", "hMuTau_SR_Events;;N", 6, 1, 7)

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

h['hMuTau_lowMET'] = ROOT.TH1F ("hMuTau_lowMET", "hMuTau_lowMET; ;N", 100, 0, 100)
h['hMuTau_lowMET_dRl'] = ROOT.TH1F ("hMuTau_lowMET_dRl", "hMuTau_lowMET_dRl;;N", 100, 0, 5)
h['hMuTau_lowMET_dRj'] = ROOT.TH1F ("hMuTau_lowMET_dRj", "hMuTau_lowMET_dRj;;N", 100, 0, 5)
h['hMuTau_lowMET_TauPt'] = ROOT.TH1F ("hMuTau_lowMET_TauPt", "hMuTau_lowMET_TauPt;;N", 500, 0, 500)
h['hMuTau_lowMET_MuonPt'] = ROOT.TH1F ("hMuTau_lowMET_MuonPt", "hMuTau_lowMET_MuonPt;;N", 500, 0, 500)
h['hMuTau_lowMET_MetPt'] = ROOT.TH1F ("hMuTau_lowMET_MetPt", "hMuTau_lowMET_MetPt;;N", 100, 0, 100)
h['hMuTau_lowMET_JetPt'] = ROOT.TH1F ("hMuTau_lowMET_JetPt", "hMuTau_lowMET_JetPt;;N", 2000, 0, 2000)
h['hMuTau_lowMET_Nj'] = ROOT.TH1F ("hMuTau_lowMET_Nj", "hMuTau_lowMET_Nj;;N", 10,0,10)
h['hMuTau_lowMET_Nbj'] = ROOT.TH1F ("hMuTau_lowMET_Nbj", "hMuTau_lowMET_Nbj;;N", 5,0,5)
h['hMuTau_lowMET_Mt'] = ROOT.TH1F ("hMuTau_lowMET_Mt", "hMuTau_lowMET_Mt;;N", 500, 0, 500)

h['hMuTau_lowMET_highMt'] = ROOT.TH1F ("hMuTau_lowMET_highMt", "hMuTau_lowMET_highMt; ;N", 100, 0, 100)
h['hMuTau_lowMET_highMt_dRl'] = ROOT.TH1F ("hMuTau_lowMET_highMt_dRl", "hMuTau_lowMET_highMt_dRl;;N", 100, 0, 5)
h['hMuTau_lowMET_highMt_dRj'] = ROOT.TH1F ("hMuTau_lowMET_highMt_dRj", "hMuTau_lowMET_highMt_dRj;;N", 100, 0, 5)
h['hMuTau_lowMET_highMt_TauPt'] = ROOT.TH1F ("hMuTau_lowMET_highMt_TauPt", "hMuTau_lowMET_highMt_TauPt;;N", 500, 0, 500)
h['hMuTau_lowMET_highMt_MuonPt'] = ROOT.TH1F ("hMuTau_lowMET_highMt_MuonPt", "hMuTau_lowMET_highMt_MuonPt;;N", 500, 0, 500)
h['hMuTau_lowMET_highMt_MetPt'] = ROOT.TH1F ("hMuTau_lowMET_highMt_MetPt", "hMuTau_lowMET_highMt_MetPt;;N", 100, 0, 100)
h['hMuTau_lowMET_highMt_JetPt'] = ROOT.TH1F ("hMuTau_lowMET_highMt_JetPt", "hMuTau_lowMET_highMt_JetPt;;N", 2000, 0, 2000)
h['hMuTau_lowMET_highMt_Nj'] = ROOT.TH1F ("hMuTau_lowMET_highMt_Nj", "hMuTau_lowMET_highMt_Nj;;N", 10,0,10)
h['hMuTau_lowMET_highMt_Nbj'] = ROOT.TH1F ("hMuTau_lowMET_highMt_Nbj", "hMuTau_lowMET_highMt_Nbj;;N", 5,0,5)
h['hMuTau_lowMET_highMt_Mt'] = ROOT.TH1F ("hMuTau_lowMET_highMt_Mt", "hMuTau_lowMET_highMt_Mt;;N", 500, 0, 500)

h['hMuTau_lowMET_lowMt'] = ROOT.TH1F ("hMuTau_lowMET_lowMt", "hMuTau_lowMET_lowMt; ;N", 100, 0, 100)
h['hMuTau_lowMET_lowMt_dRl'] = ROOT.TH1F ("hMuTau_lowMET_lowMt_dRl", "hMuTau_lowMET_lowMt_dRl;;N", 100, 0, 5)
h['hMuTau_lowMET_lowMt_dRj'] = ROOT.TH1F ("hMuTau_lowMET_lowMt_dRj", "hMuTau_lowMET_lowMt_dRj;;N", 100, 0, 5)
h['hMuTau_lowMET_lowMt_TauPt'] = ROOT.TH1F ("hMuTau_lowMET_lowMt_TauPt", "hMuTau_lowMET_lowMt_TauPt;;N", 500, 0, 500)
h['hMuTau_lowMET_lowMt_MuonPt'] = ROOT.TH1F ("hMuTau_lowMET_lowMt_MuonPt", "hMuTau_lowMET_lowMt_MuonPt;;N", 500, 0, 500)
h['hMuTau_lowMET_lowMt_MetPt'] = ROOT.TH1F ("hMuTau_lowMET_lowMt_MetPt", "hMuTau_lowMET_lowMt_MetPt;;N", 100, 0, 100)
h['hMuTau_lowMET_lowMt_JetPt'] = ROOT.TH1F ("hMuTau_lowMET_lowMt_JetPt", "hMuTau_lowMET_lowMt_JetPt;;N", 2000, 0, 2000)
h['hMuTau_lowMET_lowMt_Nj'] = ROOT.TH1F ("hMuTau_lowMET_lowMt_Nj", "hMuTau_lowMET_lowMt_Nj;;N", 10,0,10)
h['hMuTau_lowMET_lowMt_Nbj'] = ROOT.TH1F ("hMuTau_lowMET_lowMt_Nbj", "hMuTau_lowMET_lowMt_Nbj;;N", 5,0,5)
h['hMuTau_lowMET_lowMt_Mt'] = ROOT.TH1F ("hMuTau_lowMET_lowMt_Mt", "hMuTau_lowMET_lowMt_Mt;;N", 500, 0, 500)

h['hMuTau_lowMET_dRcut_MuonPt'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_MuonPt", "hMuTau_lowMET_dRcut_MuonPt;;N", 500, 0, 500)
h['hMuTau_lowMET_dRcut_TauPt'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_TauPt", "hMuTau_lowMET_dRcut_TauPt;;N", 500, 0, 500)
h['hMuTau_lowMET_dRcut'] = ROOT.TH1F ("hMuTau_lowMET_dRcut", "hMuTau_lowMET_dRcut;;N", 100, 0, 100)
h['hMuTau_lowMET_dRcut_Nj'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_Nj", "hMuTau_lowMET_dRcut_Nj;;N", 10,0,10)
h['hMuTau_lowMET_dRcut_MetPt'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_MetPt", "hMuTau_lowMET_dRcut_MetPt;;N", 100, 0, 100)
h['hMuTau_lowMET_dRcut_JetPt'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_JetPt", "hMuTau_lowMET_dRcut_JetPt;;N", 2000, 0, 2000)
h['hMuTau_lowMET_dRcut_Nbj'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_Nbj", "hMuTau_lowMET_dRcut_Nbj;;N", 5,0,5)
h['hMuTau_lowMET_dRcut_Mt'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_Mt", "hMuTau_lowMET_dRcut_Mt;;N", 500, 0, 500)

h['hMuTau_lowMET_dRcut_lowMt_TauPt'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_lowMt_TauPt", "hMuTau_lowMET_dRcut_lowMt_TauPt;;N", 500, 0, 500)
h['hMuTau_lowMET_dRcut_lowMt'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_lowMt", "hMuTau_lowMET_dRcut_lowMt;;N", 100, 0, 100)
h['hMuTau_lowMET_dRcut_lowMt_MuonPt'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_lowMt_MuonPt", "hMuTau_lowMET_dRcut_lowMt_MuonPt;;N", 500, 0, 500)
h['hMuTau_lowMET_dRcut_lowMt_MetPt'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_lowMt_MetPt", "hMuTau_lowMET_dRcut_lowMt_MetPt;;N", 100, 0, 100)
h['hMuTau_lowMET_dRcut_lowMt_JetPt'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_lowMt_JetPt", "hMuTau_lowMET_dRcut_lowMt_JetPt;;N", 2000, 0, 2000)
h['hMuTau_lowMET_dRcut_lowMt_Nj'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_lowMt_Nj", "hMuTau_lowMET_dRcut_lowMt_Nj;;N", 10,0,10)
h['hMuTau_lowMET_dRcut_lowMt_Nbj'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_lowMt_Nbj", "hMuTau_lowMET_dRcut_lowMt_Nbj;;N", 5,0,5)

h['hMuTau_lowMET_dRcut_lowMt_TauPt0'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_lowMt_TauPt0", "hMuTau_lowMET_dRcut_lowMt_TauPt0;;N", 500, 0, 500)
h['hMuTau_lowMET_dRcut_lowMt_TauPt1'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_lowMt_TauPt1", "hMuTau_lowMET_dRcut_lowMt_TauPt1;;N", 500, 0, 500)
h['hMuTau_lowMET_dRcut_lowMt_TauPt10'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_lowMt_TauPt10", "hMuTau_lowMET_dRcut_lowMt_TauPt10;;N", 500, 0, 500)

h['hMuTau_lowMET_dRcut_highMt_TauPt'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_highMt_TauPt", "hMuTau_lowMET_dRcut_highMt_TauPt;;N", 500, 0, 500)
h['hMuTau_lowMET_dRcut_highMt'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_highMt", "hMuTau_lowMET_dRcut_highMt;;N", 100, 0, 100)
h['hMuTau_lowMET_dRcut_highMt_MuonPt'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_highMt_MuonPt", "hMuTau_lowMET_dRcut_highMt_MuonPt;;N", 500, 0, 500)
h['hMuTau_lowMET_dRcut_highMt_MetPt'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_highMt_MetPt", "hMuTau_lowMET_dRcut_highMt_MetPt;;N", 100, 0, 100)
h['hMuTau_lowMET_dRcut_highMt_JetPt'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_highMt_JetPt", "hMuTau_lowMET_dRcut_highMt_JetPt;;N", 2000, 0, 2000)
h['hMuTau_lowMET_dRcut_highMt_Nj'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_highMt_Nj", "hMuTau_lowMET_dRcut_highMt_Nj;;N", 10,0,10)
h['hMuTau_lowMET_dRcut_highMt_Nbj'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_highMt_Nbj", "hMuTau_lowMET_dRcut_highMt_Nbj;;N", 5,0,5)

h['hMuTau_lowMET_dRcut_highMt_TauPt0'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_highMt_TauPt0", "hMuTau_lowMET_dRcut_highMt_TauPt0;;N", 500, 0, 500)
h['hMuTau_lowMET_dRcut_highMt_TauPt1'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_highMt_TauPt1", "hMuTau_lowMET_dRcut_highMt_TauPt1;;N", 500, 0, 500)
h['hMuTau_lowMET_dRcut_highMt_TauPt10'] = ROOT.TH1F ("hMuTau_lowMET_dRcut_highMt_TauPt10", "hMuTau_lowMET_dRcut_highMt_TauPt10;;N", 500, 0, 500)

h['hMuTau_highMt'] = ROOT.TH1F ("hMuTau_highMt", "hMuTau_highMt;;N", 100, 0, 100)
h['hMuTau_highMt_dRl'] = ROOT.TH1F ("hMuTau_highMt_dRl", "hMuTau_highMt_dRl;;N", 100, 0, 5)
h['hMuTau_highMt_dRj'] = ROOT.TH1F ("hMuTau_highMt_dRj", "hMuTau_highMt_dRj;;N", 100, 0, 5)
h['hMuTau_highMt_TauPt'] = ROOT.TH1F ("hMuTau_highMt_TauPt", "hMuTau_highMt_TauPt;;N", 500, 0, 500)
h['hMuTau_highMt_MetPt'] = ROOT.TH1F ("hMuTau_highMt_MetPt", "hMuTau_highMt_MetPt;;N", 500, 0, 500)
h['hMuTau_highMt_JetPt'] = ROOT.TH1F ("hMuTau_highMt_JetPt", "hMuTau_highMt_JetPt;;N", 2000, 0, 2000)
h['hMuTau_highMt_MuonPt'] = ROOT.TH1F ("hMuTau_highMt_MuonPt", "hMuTau_highMt_MuonPt;;N", 500, 0, 500)
h['hMuTau_highMt_Nj'] = ROOT.TH1F ("hMuTau_highMt_Nj", "hMuTau_highMt_Nj;;N", 10, 0, 10)
h['hMuTau_highMt_Nbj'] = ROOT.TH1F ("hMuTau_highMt_Nbj", "hMuTau_highMt_Nbj;;N", 5, 0, 5)

h['hMuTau_highMt_highMET'] = ROOT.TH1F ("hMuTau_highMt_highMET", "hMuTau_highMt_highMET;;N", 100, 0, 100)
h['hMuTau_highMt_highMET_dRl'] = ROOT.TH1F ("hMuTau_highMt_highMET_dRl", "hMuTau_highMt_highMET_dRl;;N", 100, 0, 5)
h['hMuTau_highMt_highMET_dRj'] = ROOT.TH1F ("hMuTau_highMt_highMET_dRj", "hMuTau_highMt_highMET_dRj;;N", 100, 0, 5)
h['hMuTau_highMt_highMET_TauPt'] = ROOT.TH1F ("hMuTau_highMt_highMET_TauPt", "hMuTau_highMt_highMET_TauPt;;N", 500, 0, 500)
h['hMuTau_highMt_highMET_MetPt'] = ROOT.TH1F ("hMuTau_highMt_highMET_MetPt", "hMuTau_highMt_highMET_MetPt;;N", 500, 0, 500)
h['hMuTau_highMt_highMET_JetPt'] = ROOT.TH1F ("hMuTau_highMt_highMET_JetPt", "hMuTau_highMt_highMET_JetPt;;N", 2000, 0, 2000)
h['hMuTau_highMt_highMET_MuonPt'] = ROOT.TH1F ("hMuTau_highMt_highMET_MuonPt", "hMuTau_highMt_highMET_MuonPt;;N", 500, 0, 500)
h['hMuTau_highMt_highMET_Nj'] = ROOT.TH1F ("hMuTau_highMt_highMET_Nj", "hMuTau_highMt_highMET_Nj;;N", 10, 0, 10)
h['hMuTau_highMt_highMET_Nbj'] = ROOT.TH1F ("hMuTau_highMt_highMET_Nbj", "hMuTau_highMt_highMET_Nbj;;N", 5, 0, 5)

h['hMuTau_highMt_dRcut_TauPt'] = ROOT.TH1F ("hMuTau_highMt_dRcut_TauPt", "hMuTau_highMt_dRcut_TauPt;;N", 500, 0, 500)
h['hMuTau_highMt_dRcut'] = ROOT.TH1F ("hMuTau_highMt_dRcut", "hMuTau_highMt_dRcut;;N", 100, 0, 100)
h['hMuTau_highMt_dRcut_Nj'] = ROOT.TH1F ("hMuTau_highMt_dRcut_Nj", "hMuTau_highMt_dRcut_Nj;;N", 10, 0, 10)
h['hMuTau_highMt_dRcut_Nbj'] = ROOT.TH1F ("hMuTau_highMt_dRcut_Nbj", "hMuTau_highMt_dRcut_Nbj;;N", 5, 0, 5)
h['hMuTau_highMt_dRcut_MetPt'] = ROOT.TH1F ("hMuTau_highMt_dRcut_MetPt", "hMuTau_highMt_dRcut_MetPt;;N", 500, 0, 500)
h['hMuTau_highMt_dRcut_JetPt'] = ROOT.TH1F ("hMuTau_highMt_dRcut_JetPt", "hMuTau_highMt_dRcut_JetPt;;N", 2000, 0, 2000)
h['hMuTau_highMt_dRcut_MuonPt'] = ROOT.TH1F ("hMuTau_highMt_dRcut_MuonPt", "hMuTau_highMt_dRcut_MuonPt;;N", 500, 0, 500)

h['hMuTau_highMt_dRcut_highMET_TauPt'] = ROOT.TH1F ("hMuTau_highMt_dRcut_highMET_TauPt", "hMuTau_highMt_dRcut_highMET_TauPt;;N", 500, 0,500)
h['hMuTau_highMt_dRcut_highMET_JetPt'] = ROOT.TH1F ("hMuTau_highMt_dRcut_highMET_JetPt", "hMuTau_highMt_dRcut_highMET_JetPt;;N", 2000, 0, 2000)
h['hMuTau_highMt_dRcut_highMET'] = ROOT.TH1F ("hMuTau_highMt_dRcut_highMET", "hMuTau_highMt_dRcut_highMET;;N", 100, 0, 100)
h['hMuTau_highMt_dRcut_highMET_Nj'] = ROOT.TH1F ("hMuTau_highMt_dRcut_highMET_Nj", "hMuTau_highMt_dRcut_highMET_Nj;;N", 10, 0, 10)
h['hMuTau_highMt_dRcut_highMET_Nbj'] = ROOT.TH1F ("hMuTau_highMt_dRcut_highMET_Nbj", "hMuTau_highMt_dRcut_highMET_Nbj;;N", 5, 0, 5)
h['hMuTau_highMt_dRcut_highMET_MuonPt'] = ROOT.TH1F ("hMuTau_highMt_dRcut_highMET_MuonPt", "hMuTau_highMt_dRcut_highMET_MuonPt;;N", 500, 0,500)

h['hMuTau_highMt_dRcut_highMET_TauPt0'] = ROOT.TH1F ("hMuTau_highMt_dRcut_highMET_TauPt0", "hMuTau_highMt_dRcut_highMET_TauPt0;;N", 500, 0, 500)
h['hMuTau_highMt_dRcut_highMET_TauPt1'] = ROOT.TH1F ("hMuTau_highMt_dRcut_highMET_TauPt1", "hMuTau_highMt_dRcut_highMET_TauPt1;;N", 500, 0, 500)
h['hMuTau_highMt_dRcut_highMET_TauPt10'] = ROOT.TH1F ("hMuTau_highMt_dRcut_highMET_TauPt10", "hMuTau_highMt_dRcut_highMET_TauPt10;;N", 500, 0, 500)

h['hMuTau_OS_TauPt'] = ROOT.TH1F ("hMuTau_OS_TauPt", "hMuTau_OS_TauPt;;N", 500, 0, 500)
h['hMuTau_OS'] = ROOT.TH1F ("hMuTau_OS", "hMuTau_OS;;N", 100, 0, 100)
h['hMuTau_OS_Nj'] = ROOT.TH1F ("hMuTau_OS_Nj", "hMuTau_OS_Nj;;N", 10, 0, 10)

h['hMuTau_SS_TauPt'] = ROOT.TH1F ("hMuTau_SS_TauPt", "hMuTau_SS_TauPt;;N", 500, 0, 500)
h['hMuTau_SS_MuonPt'] = ROOT.TH1F ("hMuTau_SS_MuonPt", "hMuTau_SS_MuonPt;;N", 500, 0, 500)
h['hMuTau_SS_MetPt'] = ROOT.TH1F ("hMuTau_SS_MetPt", "hMuTau_SS_MetPt;;N", 500, 0, 500)
h['hMuTau_SS_JetPt'] = ROOT.TH1F ("hMuTau_SS_JetPt", "hMuTau_SS_JetPt;;N", 2000, 0, 2000)
h['hMuTau_SS'] = ROOT.TH1F ("hMuTau_SS", "hMuTau_SS;;N", 100, 0, 100)
h['hMuTau_SS_Nj'] = ROOT.TH1F ("hMuTau_SS_Nj", "hMuTau_SS_Nj;;N", 10, 0, 10)
h['hMuTau_SS_Mt'] = ROOT.TH1F ("hMuTau_SS_Mt", "hMuTau_SS_Mt;;N", 500, 0, 500)
h['hMuTau_SS_dRl'] = ROOT.TH1F ("hMuTau_SS_dRl", "hMuTau_SS_dRl;;N", 100, 0, 5)
h['hMuTau_SS_dRj'] = ROOT.TH1F ("hMuTau_SS_dRj", "hMuTau_SS_dRj;;N", 100, 0, 5)

h['hMuTau_SS_lowMET_TauPt'] = ROOT.TH1F ("hMuTau_SS_lowMET_TauPt", "hMuTau_SS_lowMET_TauPt;;N", 500, 0, 500)
h['hMuTau_SS_lowMET_MuonPt'] = ROOT.TH1F ("hMuTau_SS_lowMET_MuonPt", "hMuTau_SS_lowMET_MuonPt;;N", 500, 0, 500)
h['hMuTau_SS_lowMET_MetPt'] = ROOT.TH1F ("hMuTau_SS_lowMET_MetPt", "hMuTau_SS_lowMET_MetPt;;N", 100, 0, 100)
h['hMuTau_SS_lowMET_JetPt'] = ROOT.TH1F ("hMuTau_SS_lowMET_JetPt", "hMuTau_SS_lowMET_JetPt;;N", 2000, 0, 2000)
h['hMuTau_SS_lowMET'] = ROOT.TH1F ("hMuTau_SS_lowMET", "hMuTau_SS_lowMET;;N", 100, 0, 100)
h['hMuTau_SS_lowMET_Nj'] = ROOT.TH1F ("hMuTau_SS_lowMET_Nj", "hMuTau_SS_lowMET_Nj;;N", 10, 0, 10)
h['hMuTau_SS_lowMET_Nbj'] = ROOT.TH1F ("hMuTau_SS_lowMET_Nbj", "hMuTau_SS_lowMET_Nbj;;N", 5, 0, 5)
h['hMuTau_SS_lowMET_Mt'] = ROOT.TH1F ("hMuTau_SS_lowMET_Mt", "hMuTau_SS_lowMET_Mt;;N", 500, 0, 500)
h['hMuTau_SS_lowMET_dRl'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRl", "hMuTau_SS_lowMET_dRl;;N", 100, 0, 5)
h['hMuTau_SS_lowMET_dRj'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRj", "hMuTau_SS_lowMET_dRj;;N", 100, 0, 5)

h['hMuTau_SS_lowMET_highMt_TauPt'] = ROOT.TH1F ("hMuTau_SS_lowMET_highMt_TauPt", "hMuTau_SS_lowMET_highMt_TauPt;;N", 500, 0, 500)
h['hMuTau_SS_lowMET_highMt_MuonPt'] = ROOT.TH1F ("hMuTau_SS_lowMET_highMt_MuonPt", "hMuTau_SS_lowMET_highMt_MuonPt;;N", 500, 0, 500)
h['hMuTau_SS_lowMET_highMt_MetPt'] = ROOT.TH1F ("hMuTau_SS_lowMET_highMt_MetPt", "hMuTau_SS_lowMET_highMt_MetPt;;N", 100, 0, 100)
h['hMuTau_SS_lowMET_highMt_JetPt'] = ROOT.TH1F ("hMuTau_SS_lowMET_highMt_JetPt", "hMuTau_SS_lowMET_highMt_JetPt;;N", 2000, 0, 2000)
h['hMuTau_SS_lowMET_highMt'] = ROOT.TH1F ("hMuTau_SS_lowMET_highMt", "hMuTau_SS_lowMET_highMt;;N", 100, 0, 100)
h['hMuTau_SS_lowMET_highMt_Nj'] = ROOT.TH1F ("hMuTau_SS_lowMET_highMt_Nj", "hMuTau_SS_lowMET_highMt_Nj;;N", 10, 0, 10)
h['hMuTau_SS_lowMET_highMt_Nbj'] = ROOT.TH1F ("hMuTau_SS_lowMET_highMt_Nbj", "hMuTau_SS_lowMET_highMt_Nbj;;N", 5, 0, 5)
h['hMuTau_SS_lowMET_highMt_Mt'] = ROOT.TH1F ("hMuTau_SS_lowMET_highMt_Mt", "hMuTau_SS_lowMET_highMt_Mt;;N", 500, 0, 500)
h['hMuTau_SS_lowMET_highMt_dRl'] = ROOT.TH1F ("hMuTau_SS_lowMET_highMt_dRl", "hMuTau_SS_lowMET_highMt_dRl;;N", 100, 0, 5)
h['hMuTau_SS_lowMET_highMt_dRj'] = ROOT.TH1F ("hMuTau_SS_lowMET_highMt_dRj", "hMuTau_SS_lowMET_highMt_dRj;;N", 100, 0, 5)

h['hMuTau_SS_lowMET_lowMt_TauPt'] = ROOT.TH1F ("hMuTau_SS_lowMET_lowMt_TauPt", "hMuTau_SS_lowMET_lowMt_TauPt;;N", 500, 0, 500)
h['hMuTau_SS_lowMET_lowMt_MuonPt'] = ROOT.TH1F ("hMuTau_SS_lowMET_lowMt_MuonPt", "hMuTau_SS_lowMET_lowMt_MuonPt;;N", 500, 0, 500)
h['hMuTau_SS_lowMET_lowMt_MetPt'] = ROOT.TH1F ("hMuTau_SS_lowMET_lowMt_MetPt", "hMuTau_SS_lowMET_lowMt_MetPt;;N", 100, 0, 100)
h['hMuTau_SS_lowMET_lowMt_JetPt'] = ROOT.TH1F ("hMuTau_SS_lowMET_lowMt_JetPt", "hMuTau_SS_lowMET_lowMt_JetPt;;N", 2000, 0, 2000)
h['hMuTau_SS_lowMET_lowMt'] = ROOT.TH1F ("hMuTau_SS_lowMET_lowMt", "hMuTau_SS_lowMET_lowMt;;N", 100, 0, 100)
h['hMuTau_SS_lowMET_lowMt_Nj'] = ROOT.TH1F ("hMuTau_SS_lowMET_lowMt_Nj", "hMuTau_SS_lowMET_lowMt_Nj;;N", 10, 0, 10)
h['hMuTau_SS_lowMET_lowMt_Nbj'] = ROOT.TH1F ("hMuTau_SS_lowMET_lowMt_Nbj", "hMuTau_SS_lowMET_lowMt_Nbj;;N", 5, 0, 5)
h['hMuTau_SS_lowMET_lowMt_Mt'] = ROOT.TH1F ("hMuTau_SS_lowMET_lowMt_Mt", "hMuTau_SS_lowMET_lowMt_Mt;;N", 500, 0, 500)
h['hMuTau_SS_lowMET_lowMt_dRl'] = ROOT.TH1F ("hMuTau_SS_lowMET_lowMt_dRl", "hMuTau_SS_lowMET_lowMt_dRl;;N", 100, 0, 5)
h['hMuTau_SS_lowMET_lowMt_dRj'] = ROOT.TH1F ("hMuTau_SS_lowMET_lowMt_dRj", "hMuTau_SS_lowMET_lowMt_dRj;;N", 100, 0, 5)

h['hMuTau_SS_lowMET_dRcut_TauPt'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_TauPt", "hMuTau_SS_lowMET_dRcut_TauPt;;N", 500, 0, 500)
h['hMuTau_SS_lowMET_dRcut_MuonPt'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_MuonPt", "hMuTau_SS_lowMET_dRcut_MuonPt;;N", 500, 0, 500)
h['hMuTau_SS_lowMET_dRcut_MetPt'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_MetPt", "hMuTau_SS_lowMET_dRcut_MetPt;;N", 100, 0, 100)
h['hMuTau_SS_lowMET_dRcut_JetPt'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_JetPt", "hMuTau_SS_lowMET_dRcut_JetPt;;N", 2000, 0, 2000)
h['hMuTau_SS_lowMET_dRcut'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut", "hMuTau_SS_lowMET_dRcut;;N", 100, 0, 100)
h['hMuTau_SS_lowMET_dRcut_Nj'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_Nj", "hMuTau_SS_lowMET_dRcut_Nj;;N", 10, 0, 10)
h['hMuTau_SS_lowMET_dRcut_Nbj'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_Nbj", "hMuTau_SS_lowMET_dRcut_Nbj;;N", 5, 0, 5)
h['hMuTau_SS_lowMET_dRcut_Mt'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_Mt", "hMuTau_SS_lowMET_dRcut_Mt;;N", 500, 0, 500)

h['hMuTau_SS_lowMET_dRcut_highMt_TauPt'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_highMt_TauPt", "hMuTau_SS_lowMET_dRcut_highMt_TauPt;;N", 500, 0, 500)
h['hMuTau_SS_lowMET_dRcut_highMt_MuonPt'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_highMt_MuonPt", "hMuTau_SS_lowMET_dRcut_highMt_MuonPt;;N", 500, 0, 500)
h['hMuTau_SS_lowMET_dRcut_highMt_MetPt'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_highMt_MetPt", "hMuTau_SS_lowMET_dRcut_highMt_MetPt;;N", 100, 0, 100)
h['hMuTau_SS_lowMET_dRcut_highMt_JetPt'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_highMt_JetPt", "hMuTau_SS_lowMET_dRcut_highMt_JetPt;;N", 2000, 0, 2000)
h['hMuTau_SS_lowMET_dRcut_highMt'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_highMt", "hMuTau_SS_lowMET_dRcut_highMt;;N", 100, 0, 100)
h['hMuTau_SS_lowMET_dRcut_highMt_Nj'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_highMt_Nj", "hMuTau_SS_lowMET_dRcut_highMt_Nj;;N", 10, 0, 10)
h['hMuTau_SS_lowMET_dRcut_highMt_Nbj'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_highMt_Nbj", "hMuTau_SS_lowMET_dRcut_highMt_Nbj;;N", 5, 0, 5)
h['hMuTau_SS_lowMET_dRcut_highMt_Mt'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_highMt_Mt", "hMuTau_SS_lowMET_dRcut_highMt_Mt;;N", 500, 0, 500)

h['hMuTau_SS_lowMET_dRcut_highMt_TauPt0'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_highMt_TauPt0", "hMuTau_SS_lowMET_dRcut_highMt_TauPt0;;N", 500, 0, 500)
h['hMuTau_SS_lowMET_dRcut_highMt_TauPt1'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_highMt_TauPt1", "hMuTau_SS_lowMET_dRcut_highMt_TauPt1;;N", 500, 0, 500)
h['hMuTau_SS_lowMET_dRcut_highMt_TauPt10'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_highMt_TauPt10", "hMuTau_SS_lowMET_dRcut_highMt_TauPt10;;N", 500, 0, 500)

h['hMuTau_SS_lowMET_dRcut_lowMt_TauPt'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_lowMt_TauPt", "hMuTau_SS_lowMET_dRcut_lowMt_TauPt;;N", 500, 0, 500)
h['hMuTau_SS_lowMET_dRcut_lowMt_MuonPt'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_lowMt_MuonPt", "hMuTau_SS_lowMET_dRcut_lowMt_MuonPt;;N", 500, 0, 500)
h['hMuTau_SS_lowMET_dRcut_lowMt_MetPt'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_lowMt_MetPt", "hMuTau_SS_lowMET_dRcut_lowMt_MetPt;;N", 100, 0, 100)
h['hMuTau_SS_lowMET_dRcut_lowMt_JetPt'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_lowMt_JetPt", "hMuTau_SS_lowMET_dRcut_lowMt_JetPt;;N", 2000, 0, 2000)
h['hMuTau_SS_lowMET_dRcut_lowMt'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_lowMt", "hMuTau_SS_lowMET_dRcut_lowMt;;N", 100, 0, 100)
h['hMuTau_SS_lowMET_dRcut_lowMt_Nj'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_lowMt_Nj", "hMuTau_SS_lowMET_dRcut_lowMt_Nj;;N", 10, 0, 10)
h['hMuTau_SS_lowMET_dRcut_lowMt_Nbj'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_lowMt_Nbj", "hMuTau_SS_lowMET_dRcut_lowMt_Nbj;;N", 5, 0, 5)
h['hMuTau_SS_lowMET_dRcut_lowMt_Mt'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_lowMt_Mt", "hMuTau_SS_lowMET_dRcut_lowMt_Mt;;N", 500, 0, 500)

h['hMuTau_SS_lowMET_dRcut_lowMt_TauPt0'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_lowMt_TauPt0", "hMuTau_SS_lowMET_dRcut_lowMt_TauPt0;;N", 500, 0, 500)
h['hMuTau_SS_lowMET_dRcut_lowMt_TauPt1'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_lowMt_TauPt1", "hMuTau_SS_lowMET_dRcut_lowMt_TauPt1;;N", 500, 0, 500)
h['hMuTau_SS_lowMET_dRcut_lowMt_TauPt10'] = ROOT.TH1F ("hMuTau_SS_lowMET_dRcut_lowMt_TauPt10", "hMuTau_SS_lowMET_dRcut_lowMt_TauPt10;;N", 500, 0, 500)

h['hMuTau_SS_dRcut_TauPt'] = ROOT.TH1F ("hMuTau_SS_dRcut_TauPt", "hMuTau_SS_dRcut_TauPt;;N", 500, 0, 500)
h['hMuTau_SS_dRcut'] = ROOT.TH1F ("hMuTau_SS_dRcut", "hMuTau_SS_dRcut;;N", 100, 0, 100)
h['hMuTau_SS_dRcut_Nj'] = ROOT.TH1F ("hMuTau_SS_dRcut_Nj", "hMuTau_SS_dRcut_Nj;;N", 10, 0, 10)
h['hMuTau_SS_dRcut_Nbj'] = ROOT.TH1F ("hMuTau_SS_dRcut_Nbj", "hMuTau_SS_dRcut_Nbj;;N", 5, 0, 5)
h['hMuTau_SS_dRcut_MetPt'] = ROOT.TH1F ("hMuTau_SS_dRcut_MetPt", "hMuTau_SS_dRcut_MetPt;;N", 500, 0, 500)
h['hMuTau_SS_dRcut_JetPt'] = ROOT.TH1F ("hMuTau_SS_dRcut_JetPt", "hMuTau_SS_dRcut_JetPt;;N", 2000, 0, 2000)
h['hMuTau_SS_dRcut_MuonPt'] = ROOT.TH1F ("hMuTau_SS_dRcut_MuonPt", "hMuTau_SS_dRcut_MuonPt;;N", 500, 0, 500)
h['hMuTau_SS_dRcut_Mt'] = ROOT.TH1F ("hMuTau_SS_dRcut_Mt", "hMuTau_SS_dRcut_Mt;;N", 500, 0, 500)

h['hMuTau_SS_dRcut_highMET_TauPt'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_TauPt", "hMuTau_SS_dRcut_highMET_TauPt;;N", 500, 0, 500)
h['hMuTau_SS_dRcut_highMET'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET", "hMuTau_SS_dRcut_highMET;;N", 100, 0, 100)
h['hMuTau_SS_dRcut_highMET_Nj'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_Nj", "hMuTau_SS_dRcut_highMET_Nj;;N", 10, 0, 10)
h['hMuTau_SS_dRcut_highMET_Nbj'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_Nbj", "hMuTau_SS_dRcut_highMET_Nbj;;N", 5, 0, 5)
h['hMuTau_SS_dRcut_highMET_MetPt'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_MetPt", "hMuTau_SS_dRcut_highMET_MetPt;;N", 500, 0, 500)
h['hMuTau_SS_dRcut_highMET_JetPt'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_JetPt", "hMuTau_SS_dRcut_highMET_JetPt;;N", 2000, 0, 2000)
h['hMuTau_SS_dRcut_highMET_MuonPt'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_MuonPt", "hMuTau_SS_dRcut_highMET_MuonPt;;N", 500, 0, 500)
h['hMuTau_SS_dRcut_highMET_Mt'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_Mt", "hMuTau_SS_dRcut_highMET_Mt;;N", 500, 0, 500)

h['hMuTau_SS_dRcut_highMET_highMt_TauPt'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_highMt_TauPt", "hMuTau_SS_dRcut_highMET_highMt_TauPt;;N", 500, 0, 500)
h['hMuTau_SS_dRcut_highMET_highMt'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_highMt", "hMuTau_SS_dRcut_highMET_highMt;;N", 100, 0, 100)
h['hMuTau_SS_dRcut_highMET_highMt_Nj'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_highMt_Nj", "hMuTau_SS_dRcut_highMET_highMt_Nj;;N", 10, 0, 10)
h['hMuTau_SS_dRcut_highMET_highMt_Nbj'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_highMt_Nbj", "hMuTau_SS_dRcut_highMET_highMt_Nbj;;N", 5, 0, 5)
h['hMuTau_SS_dRcut_highMET_highMt_MetPt'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_highMt_MetPt", "hMuTau_SS_dRcut_highMET_highMt_MetPt;;N", 500, 0, 500)
h['hMuTau_SS_dRcut_highMET_highMt_JetPt'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_highMt_JetPt", "hMuTau_SS_dRcut_highMET_highMt_JetPt;;N", 2000, 0, 2000)
h['hMuTau_SS_dRcut_highMET_highMt_MuonPt'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_highMt_MuonPt", "hMuTau_SS_dRcut_highMET_highMt_MuonPt;;N", 500, 0, 500)

h['hMuTau_SS_dRcut_highMET_highMt_TauPt0'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_highMt_TauPt0", "hMuTau_SS_dRcut_highMET_highMt_TauPt0;;N", 500, 0, 500)
h['hMuTau_SS_dRcut_highMET_highMt_TauPt1'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_highMt_TauPt1", "hMuTau_SS_dRcut_highMET_highMt_TauPt1;;N", 500, 0, 500)
h['hMuTau_SS_dRcut_highMET_highMt_TauPt10'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_highMt_TauPt10", "hMuTau_SS_dRcut_highMET_highMt_TauPt10;;N", 500, 0, 500)

h['hMuTau_SS_dRcut_highMET_lowMt_TauPt'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_lowMt_TauPt", "hMuTau_SS_dRcut_highMET_lowMt_TauPt;;N", 500, 0, 500)
h['hMuTau_SS_dRcut_highMET_lowMt'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_lowMt", "hMuTau_SS_dRcut_highMET_lowMt;;N", 100, 0, 100)
h['hMuTau_SS_dRcut_highMET_lowMt_Nj'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_lowMt_Nj", "hMuTau_SS_dRcut_highMET_lowMt_Nj;;N", 10, 0, 10)
h['hMuTau_SS_dRcut_highMET_lowMt_Nbj'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_lowMt_Nbj", "hMuTau_SS_dRcut_highMET_lowMt_Nbj;;N", 5, 0, 5)
h['hMuTau_SS_dRcut_highMET_lowMt_MetPt'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_lowMt_MetPt", "hMuTau_SS_dRcut_highMET_lowMt_MetPt;;N", 500, 0, 500)
h['hMuTau_SS_dRcut_highMET_lowMt_JetPt'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_lowMt_JetPt", "hMuTau_SS_dRcut_highMET_lowMt_JetPt;;N", 2000, 0, 2000)
h['hMuTau_SS_dRcut_highMET_lowMt_MuonPt'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_lowMt_MuonPt", "hMuTau_SS_dRcut_highMET_lowMt_MuonPt;;N", 500, 0, 500)

h['hMuTau_SS_dRcut_highMET_lowMt_TauPt0'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_lowMt_TauPt0", "hMuTau_SS_dRcut_highMET_lowMt_TauPt0;;N", 500, 0, 500)
h['hMuTau_SS_dRcut_highMET_lowMt_TauPt1'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_lowMt_TauPt1", "hMuTau_SS_dRcut_highMET_lowMt_TauPt1;;N", 500, 0, 500)
h['hMuTau_SS_dRcut_highMET_lowMt_TauPt10'] = ROOT.TH1F ("hMuTau_SS_dRcut_highMET_lowMt_TauPt10", "hMuTau_SS_dRcut_highMET_lowMt_TauPt10;;N", 500, 0, 500)

if isData == 0 or isAltered == 1:
    h['hMuTau_SR'] = ROOT.TH1F ("hMuTau_SR", "hMuTau_SR; ;N", 100, 0, 100)
    h['hMuTau_SR_dRl'] = ROOT.TH1F ("hMuTau_SR_dRl", "hMuTau_SR_dRl;;N", 100, 0, 5)
    h['hMuTau_SR_dRj'] = ROOT.TH1F ("hMuTau_SR_dRj", "hMuTau_SR_dRj;;N", 100, 0, 5)
    h['hMuTau_SR_TauPt'] = ROOT.TH1F ("hMuTau_SR_TauPt", "hMuTau_SR_TauPt;;N", 500, 0, 500)
    h['hMuTau_SR_MuonPt'] = ROOT.TH1F ("hMuTau_SR_MuonPt", "hMuTau_SR_MuonPt;;N", 500, 0, 500)
    h['hMuTau_SR_MetPt'] = ROOT.TH1F ("hMuTau_SR_MetPt", "hMuTau_SR_MetPt;;N", 500, 0, 500)
    h['hMuTau_SR_JetPt'] = ROOT.TH1F ("hMuTau_SR_JetPt", "hMuTau_SR_JetPt;;N", 2000, 0, 2000)
    h['hMuTau_SR_Nj'] = ROOT.TH1F ("hMuTau_SR_Nj", "hMuTau_SR_Nj;;N", 10,0,10)
    h['hMuTau_SR_Nbj'] = ROOT.TH1F ("hMuTau_SR_Nbj", "hMuTau_SR_Nbj;;N", 5,0,5)
    h['hMuTau_SR_Mt'] = ROOT.TH1F ("hMuTau_SR_Mt", "hMuTau_SR_Mt;;N", 500, 0, 500)

    h['hMuTau_SR_dRcut'] = ROOT.TH1F ("hMuTau_SR_dRcut", "hMuTau_SR_dRcut; ;N", 100, 0, 100)
    h['hMuTau_SR_dRcut_dRl'] = ROOT.TH1F ("hMuTau_SR_dRcut_dRl", "hMuTau_SR_dRcut_dRl;;N", 100, 0, 5)
    h['hMuTau_SR_dRcut_dRj'] = ROOT.TH1F ("hMuTau_SR_dRcut_dRj", "hMuTau_SR_dRcut_dRj;;N", 100, 0, 5)
    h['hMuTau_SR_dRcut_TauPt'] = ROOT.TH1F ("hMuTau_SR_dRcut_TauPt", "hMuTau_SR_dRcut_TauPt;;N", 500, 0, 500)
    h['hMuTau_SR_dRcut_MuonPt'] = ROOT.TH1F ("hMuTau_SR_dRcut_MuonPt", "hMuTau_SR_dRcut_MuonPt;;N", 500, 0, 500)
    h['hMuTau_SR_dRcut_MetPt'] = ROOT.TH1F ("hMuTau_SR_dRcut_MetPt", "hMuTau_SR_dRcut_MetPt;;N", 500, 0, 500)
    h['hMuTau_SR_dRcut_JetPt'] = ROOT.TH1F ("hMuTau_SR_dRcut_JetPt", "hMuTau_SR_dRcut_JetPt;;N", 2000, 0, 2000)
    h['hMuTau_SR_dRcut_Nj'] = ROOT.TH1F ("hMuTau_SR_dRcut_Nj", "hMuTau_SR_dRcut_Nj;;N", 10,0,10)
    h['hMuTau_SR_dRcut_Nbj'] = ROOT.TH1F ("hMuTau_SR_dRcut_Nbj", "hMuTau_SR_dRcut_Nbj;;N", 5,0,5)
    h['hMuTau_SR_dRcut_Mt'] = ROOT.TH1F ("hMuTau_SR_dRcut_Mt", "hMuTau_SR_dRcut_Mt;;N", 500, 0, 500)

    h['hMuTau_SR_dRcut_highMET'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET", "hMuTau_SR_dRcut_highMET; ;N", 100, 0, 100)
    h['hMuTau_SR_dRcut_highMET_dRl'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_dRl", "hMuTau_SR_dRcut_highMET_dRl;;N", 100, 0, 5)
    h['hMuTau_SR_dRcut_highMET_dRj'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_dRj", "hMuTau_SR_dRcut_highMET_dRj;;N", 100, 0, 5)
    h['hMuTau_SR_dRcut_highMET_TauPt'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_TauPt", "hMuTau_SR_dRcut_highMET_TauPt;;N", 500, 0, 500)
    h['hMuTau_SR_dRcut_highMET_MuonPt'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_MuonPt", "hMuTau_SR_dRcut_highMET_MuonPt;;N", 500, 0, 500)
    h['hMuTau_SR_dRcut_highMET_MetPt'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_MetPt", "hMuTau_SR_dRcut_highMET_MetPt;;N", 500, 0, 500)
    h['hMuTau_SR_dRcut_highMET_JetPt'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_JetPt", "hMuTau_SR_dRcut_highMET_JetPt;;N", 2000, 0, 2000)
    h['hMuTau_SR_dRcut_highMET_Nj'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_Nj", "hMuTau_SR_dRcut_highMET_Nj;;N", 10,0,10)
    h['hMuTau_SR_dRcut_highMET_Nbj'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_Nbj", "hMuTau_SR_dRcut_highMET_Nbj;;N", 5,0,5)
    h['hMuTau_SR_dRcut_highMET_Mt'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_Mt", "hMuTau_SR_dRcut_highMET_Mt;;N", 500, 0, 500)

    h['hMuTau_SR_dRcut_highMET_lowMt'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_lowMt", "hMuTau_SR_dRcut_highMET_lowMt; ;N", 100, 0, 100)
    h['hMuTau_SR_dRcut_highMET_lowMt_dRl'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_lowMt_dRl", "hMuTau_SR_dRcut_highMET_lowMt_dRl;;N", 100, 0, 5)
    h['hMuTau_SR_dRcut_highMET_lowMt_dRj'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_lowMt_dRj", "hMuTau_SR_dRcut_highMET_lowMt_dRj;;N", 100, 0, 5)
    h['hMuTau_SR_dRcut_highMET_lowMt_TauPt'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_lowMt_TauPt", "hMuTau_SR_dRcut_highMET_lowMt_TauPt;;N", 500, 0, 500)
    h['hMuTau_SR_dRcut_highMET_lowMt_MuonPt'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_lowMt_MuonPt", "hMuTau_SR_dRcut_highMET_lowMt_MuonPt;;N", 500, 0, 500)
    h['hMuTau_SR_dRcut_highMET_lowMt_MetPt'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_lowMt_MetPt", "hMuTau_SR_dRcut_highMET_lowMt_MetPt;;N", 500, 0, 500)
    h['hMuTau_SR_dRcut_highMET_lowMt_JetPt'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_lowMt_JetPt", "hMuTau_SR_dRcut_highMET_lowMt_JetPt;;N", 2000, 0, 2000)
    h['hMuTau_SR_dRcut_highMET_lowMt_Nj'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_lowMt_Nj", "hMuTau_SR_dRcut_highMET_lowMt_Nj;;N", 10,0,10)
    h['hMuTau_SR_dRcut_highMET_lowMt_Nbj'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_lowMt_Nbj", "hMuTau_SR_dRcut_highMET_lowMt_Nbj;;N", 5,0,5)
    h['hMuTau_SR_dRcut_highMET_lowMt_Mt'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_lowMt_Mt", "hMuTau_SR_dRcut_highMET_lowMt_Mt;;N", 500, 0, 500)

    h['hMuTau_SR_dRcut_highMET_lowMt_TauPt0'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_lowMt_TauPt0", "hMuTau_SR_dRcut_highMET_lowMt_TauPt0;;N", 500, 0, 500)
    h['hMuTau_SR_dRcut_highMET_lowMt_TauPt1'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_lowMt_TauPt1", "hMuTau_SR_dRcut_highMET_lowMt_TauPt1;;N", 500, 0, 500)
    h['hMuTau_SR_dRcut_highMET_lowMt_TauPt10'] = ROOT.TH1F ("hMuTau_SR_dRcut_highMET_lowMt_TauPt10", "hMuTau_SR_dRcut_highMET_lowMt_TauPt10;;N", 500, 0, 500)


if plot2DforTau == 1:
    h['hMuTau_SS_dRcut_highMET_lowMt_TauPtJetPt'] = ROOT.TH2F("hMuTau_SS_dRcut_highMET_lowMt_TauPtJetPt","hMuTau_SS_dRcut_highMET_lowMt_TauPtJetPt;TauPt;JetPt",500,0,500,2000,0,2000)
    h['hMuTau_SS_dRcut_highMET_lowMt_TauPtMuonPt'] = ROOT.TH2F("hMuTau_SS_dRcut_highMET_lowMt_TauPtMuonPt","hMuTau_SS_dRcut_highMET_lowMt_TauPtMuonPt;TauPt;MuonPt",500,0,500,500,0,500)
    h['hMuTau_SS_dRcut_highMET_lowMt_TauPtMetPt'] = ROOT.TH2F("hMuTau_SS_dRcut_highMET_lowMt_TauPtMetPt","hMuTau_SS_dRcut_highMET_lowMt_TauPtMetPt;TauPt;MetPt",500,0,500,500,0,500)
    h['hMuTau_SS_dRcut_highMET_lowMt_TauPtMuMuMass'] = ROOT.TH2F("hMuTau_SS_dRcut_highMET_lowMt_TauPtMuMuMass","hMuTau_SS_dRcut_highMET_lowMt_TauPtMuMuMass;TauPt;M_{#mu#mu}",500,0,500,150,0,150)
    h['hMuTau_SS_dRcut_highMET_lowMt_TauPtdRj'] = ROOT.TH2F("hMuTau_SS_dRcut_highMET_lowMt_TauPtdRj","hMuTau_SS_dRcut_highMET_lowMt_TauPtdRj;TauPt;dR(tau,jet1)", 500, 0, 500,100,0,5)
    h['hMuTau_SS_dRcut_highMET_lowMt_TauPtdRj2'] = ROOT.TH2F("hMuTau_SS_dRcut_highMET_lowMt_TauPtdRj2","hMuTau_SS_dRcut_highMET_lowMt_TauPtdRj2;TauPt;dR(tau,jet2)",500,0,500,100,0,5)

    h['hMuTau_SS_dRcut_highMET_highMt_TauPtJetPt'] = ROOT.TH2F("hMuTau_SS_dRcut_highMET_highMt_TauPtJetPt","hMuTau_SS_dRcut_highMET_highMt_TauPtJetPt;TauPt;JetPt",500,0,500,2000,0,2000)
    h['hMuTau_SS_dRcut_highMET_highMt_TauPtMuonPt'] = ROOT.TH2F("hMuTau_SS_dRcut_highMET_highMt_TauPtMuonPt","hMuTau_SS_dRcut_highMET_highMt_TauPtMuonPt;TauPt;MuonPt",500,0,500,500,0,500)
    h['hMuTau_SS_dRcut_highMET_highMt_TauPtMetPt'] = ROOT.TH2F("hMuTau_SS_dRcut_highMET_highMt_TauPtMetPt","hMuTau_SS_dRcut_highMET_highMt_TauPtMetPt;TauPt;MetPt",500,0,500,500,0,500)
    h['hMuTau_SS_dRcut_highMET_highMt_TauPtMuMuMass'] = ROOT.TH2F("hMuTau_SS_dRcut_highMET_highMt_TauPtMuMuMass","hMuTau_SS_dRcut_highMET_highMt_TauPtMuMuMass;TauPt;M_{#mu#mu}",500,0,500,150,0,150)
    h['hMuTau_SS_dRcut_highMET_highMt_TauPtdRj'] = ROOT.TH2F("hMuTau_SS_dRcut_highMET_highMt_TauPtdRj","hMuTau_SS_dRcut_highMET_highMt_TauPtdRj;TauPt;dR(tau,jet1)", 500, 0, 500,100,0,5)
    h['hMuTau_SS_dRcut_highMET_highMt_TauPtdRj2'] = ROOT.TH2F("hMuTau_SS_dRcut_highMET_highMt_TauPtdRj2","hMuTau_SS_dRcut_highMET_highMt_TauPtdRj2;TauPt;dR(tau,jet2)",500,0,500,100,0,5)

    h['hMuTau_SS_lowMET_dRcut_highMt_TauPtJetPt'] = ROOT.TH2F("hMuTau_SS_lowMET_dRcut_highMt_TauPtJetPt","hMuTau_SS_lowMET_dRcut_highMt_TauPtJetPt;TauPt;JetPt",500,0,500,2000,0,2000)
    h['hMuTau_SS_lowMET_dRcut_highMt_TauPtMuonPt'] = ROOT.TH2F("hMuTau_SS_lowMET_dRcut_highMt_TauPtMuonPt","hMuTau_SS_lowMET_dRcut_highMt_TauPtMuonPt;TauPt;MuonPt",500,0,500,500,0,500)
    h['hMuTau_SS_lowMET_dRcut_highMt_TauPtMetPt'] = ROOT.TH2F("hMuTau_SS_lowMET_dRcut_highMt_TauPtMetPt","hMuTau_SS_lowMET_dRcut_highMt_TauPtMetPt;TauPt;MetPt",500,0,500,500,0,500)
    h['hMuTau_SS_lowMET_dRcut_highMt_TauPtMuMuMass'] = ROOT.TH2F("hMuTau_SS_lowMET_dRcut_highMt_TauPtMuMuMass","hMuTau_SS_lowMET_dRcut_highMt_TauPtMuMuMass;TauPt;M_{#mu#mu}",500,0,500,150,0,150)
    h['hMuTau_SS_lowMET_dRcut_highMt_TauPtdRj'] = ROOT.TH2F("hMuTau_SS_lowMET_dRcut_highMt_TauPtdRj","hMuTau_SS_lowMET_dRcut_highMt_TauPtdRj;TauPt;dR(tau,jet1)", 500, 0, 500,100,0,5)
    h['hMuTau_SS_lowMET_dRcut_highMt_TauPtdRj2'] = ROOT.TH2F("hMuTau_SS_lowMET_dRcut_highMt_TauPtdRj2","hMuTau_SS_lowMET_dRcut_highMt_TauPtdRj2;TauPt;dR(tau,jet2)",500,0,500,100,0,5)

    h['hMuTau_SS_lowMET_dRcut_lowMt_TauPtJetPt'] = ROOT.TH2F("hMuTau_SS_lowMET_dRcut_lowMt_TauPtJetPt","hMuTau_SS_lowMET_dRcut_lowMt_TauPtJetPt;TauPt;JetPt",500,0,500,2000,0,2000)
    h['hMuTau_SS_lowMET_dRcut_lowMt_TauPtMuonPt'] = ROOT.TH2F("hMuTau_SS_lowMET_dRcut_lowMt_TauPtMuonPt","hMuTau_SS_lowMET_dRcut_lowMt_TauPtMuonPt;TauPt;MuonPt",500,0,500,500,0,500)
    h['hMuTau_SS_lowMET_dRcut_lowMt_TauPtMetPt'] = ROOT.TH2F("hMuTau_SS_lowMET_dRcut_lowMt_TauPtMetPt","hMuTau_SS_lowMET_dRcut_lowMt_TauPtMetPt;TauPt;MetPt",500,0,500,500,0,500)
    h['hMuTau_SS_lowMET_dRcut_lowMt_TauPtMuMuMass'] = ROOT.TH2F("hMuTau_SS_lowMET_dRcut_lowMt_TauPtMuMuMass","hMuTau_SS_lowMET_dRcut_lowMt_TauPtMuMuMass;TauPt;M_{#mu#mu}",500,0,500,150,0,150)
    h['hMuTau_SS_lowMET_dRcut_lowMt_TauPtdRj'] = ROOT.TH2F("hMuTau_SS_lowMET_dRcut_lowMt_TauPtdRj","hMuTau_SS_lowMET_dRcut_lowMt_TauPtdRj;TauPt;dR(tau,jet1)", 500, 0, 500,100,0,5)
    h['hMuTau_SS_lowMET_dRcut_lowMt_TauPtdRj2'] = ROOT.TH2F("hMuTau_SS_lowMET_dRcut_lowMt_TauPtdRj2","hMuTau_SS_lowMET_dRcut_lowMt_TauPtdRj2;TauPt;dR(tau,jet2)",500,0,500,100,0,5)




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
            if met_pt > 100:
                if abs(m.DeltaPhi(e1)) < 1.0 and abs(m.DeltaPhi(j)) > 2.0:
                    if ( j.Pt() > 500 and ( isSingleJet == 1 or isHT == 1 ) ):
                        isEE = 1
                        h['hTauEvents'].Fill(1,1)
                    # if ( j.Pt() > 500 and ( isSingleJet == 1 or isHT == 1 ) ) \
                    #    or ( j.Pt() > 500 and met_pt > 200 and isHTMHT == 1 ) \
                    #    or ( e1.Pt() > 35 and isIsoEle == 1 ):
                    #     isEE = 1
                    #     h['hTauEvents'].Fill(1,1)

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

           if plotSemiLeptonic == 1:
               h['hMuMu_Baseline_Nbj'].Fill(len(s_b), genweight)
               if mu1.DeltaR(mu2) < 0.4 and mu1.DeltaR(j)> 0.8  and mu2.DeltaR(j) > 0.8:
                   h['hMuMu_dRcut_Nbj'].Fill(len(s_b), genweight)
                   if met_pt > 100:
                       if abs(m.DeltaPhi(mu1)) < 1 and abs(m.DeltaPhi(j)) > 2:
                           h['hMuMu_highMET_Nbj'].Fill(len(s_b), genweight)

           if plotSemiLeptonic == 1:
           
               if met_pt < 100 and len(s_b) == 0:
                   h['hMuMu_Events'].Fill(3, genweight)
                   if isJetHTEvent == 1 : h['hMuMu_Trigger'].Fill(7, genweight)
                   h['hMuMu_Inclusive_Trigger'].Fill(5, genweight)
                   h['hMuMu_lowMET'].Fill((mu1+mu2).M(), genweight)
                   h['hMuMu_lowMET_MetPt'].Fill(met_pt, genweight)
                   h['hMuMu_lowMET_dRl'].Fill(mu1.DeltaR(mu2), genweight)
                   h['hMuMu_lowMET_dRj'].Fill(mu1.DeltaR(j), genweight)
                   h['hMuMu_lowMET_Nj'].Fill(len(s_j), genweight)
                   h['hMuMu_lowMET_Muon1Pt'].Fill(mu1.Pt(), genweight)
                   h['hMuMu_lowMET_Muon2Pt'].Fill(mu2.Pt(), genweight)
                   h['hMuMu_lowMET_JetPt'].Fill(j.Pt(), genweight)
                   h['hMuMu_lowMET_dPhil'].Fill(m.DeltaPhi(mu1), genweight)
                   h['hMuMu_lowMET_dPhij'].Fill(m.DeltaPhi(j), genweight)

                   if mu1.DeltaR(mu2) < 0.4 and mu1.DeltaR(j)> 0.8 and mu2.DeltaR(j) > 0.8:
                       h['hMuMu_Events'].Fill(4, genweight)
                       h['hMuMu_dRcut'].Fill((mu1+mu2).M(), genweight)
                       h['hMuMu_dRcut_MetPt'].Fill(met_pt, genweight)
                       h['hMuMu_dRcut_dRl'].Fill(mu1.DeltaR(mu2), genweight)
                       h['hMuMu_dRcut_dRj'].Fill(mu1.DeltaR(j), genweight)
                       h['hMuMu_dRcut_Nj'].Fill(len(s_j), genweight)
                       h['hMuMu_lowMET_dRcut_Nbj'].Fill(len(s_b), genweight)
                       h['hMuMu_dRcut_Muon1Pt'].Fill(mu1.Pt(), genweight)
                       h['hMuMu_dRcut_Muon2Pt'].Fill(mu2.Pt(), genweight)
                       h['hMuMu_dRcut_JetPt'].Fill(j.Pt(), genweight)
                       h['hMuMu_dRcut_dPhil'].Fill(m.DeltaPhi(mu1), genweight)
                       h['hMuMu_dRcut_dPhij'].Fill(m.DeltaPhi(j), genweight)

                       if abs(m.DeltaPhi(mu1)) < 1 and abs(m.DeltaPhi(j)) > 2:
                           h['hMuMu_Events'].Fill(5, genweight)
                           h['hMuMu_dPhicut'].Fill((mu1+mu2).M(), genweight)
                           h['hMuMu_dPhicut_MetPt'].Fill(met_pt, genweight)
                           h['hMuMu_dPhicut_dRl'].Fill(mu1.DeltaR(mu2), genweight)
                           h['hMuMu_dPhicut_dRj'].Fill(mu1.DeltaR(j), genweight)
                           h['hMuMu_dPhicut_Nj'].Fill(len(s_j), genweight)
                           h['hMuMu_dPhicut_Nbj'].Fill(len(s_b), genweight)
                           h['hMuMu_dPhicut_Muon1Pt'].Fill(mu1.Pt(), genweight)
                           h['hMuMu_dPhicut_Muon2Pt'].Fill(mu2.Pt(), genweight)
                           h['hMuMu_dPhicut_JetPt'].Fill(j.Pt(), genweight)
                           h['hMuMu_dPhicut_dPhil'].Fill(m.DeltaPhi(mu1), genweight)
                           h['hMuMu_dPhicut_dPhij'].Fill(m.DeltaPhi(j), genweight)

           if len(s_b) == 0:
                if mu1.DeltaR(mu2) < 0.4 and mu1.DeltaR(j)> 0.8  and mu2.DeltaR(j) > 0.8:
                   if met_pt > 100:
                       if abs(m.DeltaPhi(mu1)) < 1 and abs(m.DeltaPhi(j)) > 2:
                           isMuMu = 1
                           h['hTauEvents'].Fill(2,genweight)
                           h['hMuMu_SR_dRcut_highMET_dPhicut'].Fill((mu1+mu2).M(), genweight)


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

   if ( isData == 0 and ( ( plotHTTrig == 1 and isJetHTEvent == 1 ) or ( plotJetHT == 1 and isJetHTEvent == 1 ) or ( plotSingleMuon == 1 and isSingleMuonEvent == 1 ) ) ) or \
      ( isData == 1 and ( isJetHTDataset == 1 and ( isJetHTEvent == 1 and isSingleMuonEvent == 0 ) ) ) or \
      ( isData == 1 and ( isSingleMuonDataset == 1 and isSingleMuonEvent == 1 ) ):

       h['hEMu_Events'].Fill(2, genweight)
       if plotSemiLeptonic == 1:

           if j.Pt() > 100 and (e+mu).M() > 1:
               if met_pt < 100:
                   h['hEMu_Events'].Fill(3, genweight)
                   if isJetHTEvent == 1 : h['hEMu_Trigger'].Fill(7, genweight)
                   h['hEMu_Inclusive_Trigger'].Fill(5, genweight)
                   h['hEMu_lowMET'].Fill((mu+e).M(), genweight)
                   h['hEMu_lowMET_MetPt'].Fill(met_pt, genweight)
                   h['hEMu_lowMET_dRl'].Fill(mu.DeltaR(e), genweight)
                   h['hEMu_lowMET_dRj'].Fill(mu.DeltaR(j), genweight)
                   h['hEMu_lowMET_Nj'].Fill(len(s_j), genweight)
                   h['hEMu_lowMET_MuonPt'].Fill(mu.Pt(), genweight)
                   h['hEMu_lowMET_ElectronPt'].Fill(e.Pt(), genweight)
                   h['hEMu_lowMET_JetPt'].Fill(j.Pt(), genweight)
                   h['hEMu_lowMET_dPhil'].Fill(m.DeltaPhi(mu), genweight)
                   h['hEMu_lowMET_dPhij'].Fill(m.DeltaPhi(j), genweight)

                   if mu.DeltaR(e) < 0.4 and mu.DeltaR(j)> 0.8  and e.DeltaR(j) > 0.8:
                       h['hEMu_Events'].Fill(4, genweight)
                       h['hEMu_dRcut'].Fill((mu+e).M(), genweight)
                       h['hEMu_dRcut_MetPt'].Fill(met_pt, genweight)
                       h['hEMu_dRcut_dRl'].Fill(mu.DeltaR(e), genweight)
                       h['hEMu_dRcut_dRj'].Fill(mu.DeltaR(j), genweight)
                       h['hEMu_dRcut_Nj'].Fill(len(s_j), genweight)
                       h['hEMu_lowMET_dRcut_Nbj'].Fill(len(s_b), genweight)
                       h['hEMu_dRcut_MuonPt'].Fill(mu.Pt(), genweight)
                       h['hEMu_dRcut_ElectronPt'].Fill(e.Pt(), genweight)
                       h['hEMu_dRcut_JetPt'].Fill(j.Pt(), genweight)
                       h['hEMu_dRcut_dPhil'].Fill(m.DeltaPhi(mu), genweight)
                       h['hEMu_dRcut_dPhij'].Fill(m.DeltaPhi(j), genweight)

       if (e+mu).M() > 1 :
           if mu.DeltaR(e) < 0.4 and mu.DeltaR(j)> 0.8  and e.DeltaR(j) > 0.8:
               if met_pt > 100.0:
                   if abs(m.DeltaPhi(mu)) < 1 and abs( m.DeltaPhi(j) ) > 2:
                       if j.Pt() > 100 :
                           isEMu = 1
                           h['hTauEvents'].Fill(3,genweight)
                           h['hEMu_SR_dRcut_highMET_dPhicut'].Fill((mu+e).M(), genweight)
                    
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

           if j.Pt() > 100 and (mu+tau).M() > 1:
                h['hMuTau_Events'].Fill(2, genweight)

                if lmuon[0].charge*mtau[0].charge < 0 :
                    h['hMuTau_OS_Nbj'].Fill(len(s_b), genweight)

                if lmuon[0].charge*mtau[0].charge < 0 and len(s_b) == 0: #OS
                    h['hMuTau_OS_TauPt'].Fill(tau.Pt(), genweight)
                    h['hMuTau_OS'].Fill((mu+tau).M(), genweight)
                    h['hMuTau_OS_Nj'].Fill(len(s_j), genweight)
                    # plot met pt

                    if met_pt < 100 : #low MET
                        h['hMuTau_Events'].Fill(3, genweight)
                        if isJetHTEvent == 1 : h['hMuTau_Trigger'].Fill(7, genweight)
                        h['hMuTau_Inclusive_Trigger'].Fill(5, genweight)
                        h['hMuTau_lowMET_dRl'].Fill(mu.DeltaR(tau), genweight)
                        h['hMuTau_lowMET_dRj'].Fill(j.DeltaR(mu+tau), genweight)
                        h['hMuTau_lowMET_TauPt'].Fill(tau.Pt(), genweight)
                        h['hMuTau_lowMET_MuonPt'].Fill(mu.Pt(), genweight)
                        h['hMuTau_lowMET_MetPt'].Fill(met_pt, genweight)
                        h['hMuTau_lowMET_JetPt'].Fill(j.Pt(), genweight)
                        h['hMuTau_lowMET'].Fill((mu+tau).M(), genweight)
                        h['hMuTau_lowMET_Nj'].Fill(len(s_j), genweight)
                        h['hMuTau_lowMET_Nbj'].Fill(len(s_b), genweight)
                        h['hMuTau_lowMET_Mt'].Fill(Mt(mu,m), genweight)
                        
                        if mu.DeltaR(tau) < 0.4 and mu.DeltaR(j) > 0.8 and tau.DeltaR(j) > 0.8 : #dR cut

                            h['hMuTau_Events'].Fill(4, genweight)
                            h['hMuTau_lowMET_dRcut_TauPt'].Fill(tau.Pt(), genweight)
                            h['hMuTau_lowMET_dRcut_MuonPt'].Fill(mu.Pt(), genweight)
                            h['hMuTau_lowMET_dRcut_MetPt'].Fill(met_pt, genweight)
                            h['hMuTau_lowMET_dRcut_JetPt'].Fill(j.Pt(), genweight)
                            h['hMuTau_lowMET_dRcut'].Fill((mu+tau).M(), genweight)
                            h['hMuTau_lowMET_dRcut_Nj'].Fill(len(s_j), genweight)
                            h['hMuTau_lowMET_dRcut_Nbj'].Fill(len(s_b), genweight)
                            h['hMuTau_lowMET_dRcut_Mt'].Fill(Mt(mu,m), genweight)

                            if Mt(mu,m) < 50 : #lowMET, dRcut, lowMt
                                h['hMuTau_Events'].Fill(5, genweight)
                                h['hMuTau_lowMET_dRcut_lowMt_TauPt'].Fill(tau.Pt(), genweight)
                                h['hMuTau_lowMET_dRcut_lowMt_MuonPt'].Fill(mu.Pt(), genweight)
                                h['hMuTau_lowMET_dRcut_lowMt_MetPt'].Fill(met_pt, genweight)
                                h['hMuTau_lowMET_dRcut_lowMt_JetPt'].Fill(j.Pt(), genweight)
                                h['hMuTau_lowMET_dRcut_lowMt'].Fill((mu+tau).M(), genweight)
                                h['hMuTau_lowMET_dRcut_lowMt_Nj'].Fill(len(s_j), genweight)
                                h['hMuTau_lowMET_dRcut_lowMt_Nbj'].Fill(len(s_b), genweight)

                                if mtau[0].decaymode == 0 : h['hMuTau_lowMET_dRcut_lowMt_TauPt0'].Fill(tau.Pt(), genweight)
                                if mtau[0].decaymode == 1 : h['hMuTau_lowMET_dRcut_lowMt_TauPt1'].Fill(tau.Pt(), genweight)
                                if mtau[0].decaymode == 10 : h['hMuTau_lowMET_dRcut_lowMt_TauPt10'].Fill(tau.Pt(), genweight)

                            if Mt(mu,m) > 50 : #lowMET, dRcut, highMt
                                h['hMuTau_lowMET_dRcut_highMt_TauPt'].Fill(tau.Pt(), genweight)
                                h['hMuTau_lowMET_dRcut_highMt_MuonPt'].Fill(mu.Pt(), genweight)
                                h['hMuTau_lowMET_dRcut_highMt_MetPt'].Fill(met_pt, genweight)
                                h['hMuTau_lowMET_dRcut_highMt_JetPt'].Fill(j.Pt(), genweight)
                                h['hMuTau_lowMET_dRcut_highMt'].Fill((mu+tau).M(), genweight)
                                h['hMuTau_lowMET_dRcut_highMt_Nj'].Fill(len(s_j), genweight)
                                h['hMuTau_lowMET_dRcut_highMt_Nbj'].Fill(len(s_b), genweight)

                                if mtau[0].decaymode == 0 : h['hMuTau_lowMET_dRcut_highMt_TauPt0'].Fill(tau.Pt(), genweight)
                                if mtau[0].decaymode == 1 : h['hMuTau_lowMET_dRcut_highMt_TauPt1'].Fill(tau.Pt(), genweight)
                                if mtau[0].decaymode == 10 : h['hMuTau_lowMET_dRcut_highMt_TauPt10'].Fill(tau.Pt(), genweight)

                        if Mt(mu,m) > 50 : #lowMET, highmT
                            h['hMuTau_lowMET_highMt_TauPt'].Fill(tau.Pt(), genweight)
                            h['hMuTau_lowMET_highMt_MuonPt'].Fill(mu.Pt(), genweight)
                            h['hMuTau_lowMET_highMt_MetPt'].Fill(met_pt, genweight)
                            h['hMuTau_lowMET_highMt_JetPt'].Fill(j.Pt(), genweight)
                            h['hMuTau_lowMET_highMt'].Fill((mu+tau).M(), genweight)
                            h['hMuTau_lowMET_highMt_Nj'].Fill(len(s_j), genweight)
                            h['hMuTau_lowMET_highMt_dRl'].Fill(mu.DeltaR(tau), genweight)
                            h['hMuTau_lowMET_highMt_dRj'].Fill(j.DeltaR(mu+tau), genweight)


                        if Mt(mu,m) < 50 : #lowMET, lowmT
                            h['hMuTau_lowMET_lowMt_TauPt'].Fill(tau.Pt(), genweight)
                            h['hMuTau_lowMET_lowMt_MuonPt'].Fill(mu.Pt(), genweight)
                            h['hMuTau_lowMET_lowMt_MetPt'].Fill(met_pt, genweight)
                            h['hMuTau_lowMET_lowMt_JetPt'].Fill(j.Pt(), genweight)
                            h['hMuTau_lowMET_lowMt'].Fill((mu+tau).M(), genweight)
                            h['hMuTau_lowMET_lowMt_Nj'].Fill(len(s_j), genweight)
                            h['hMuTau_lowMET_lowMt_dRl'].Fill(mu.DeltaR(tau), genweight)
                            h['hMuTau_lowMET_lowMt_dRj'].Fill(j.DeltaR(mu+tau), genweight)

                    if Mt(mu,m) > 50 : #high mT
                        h['hMuTau_highMt_dRl'].Fill(mu.DeltaR(tau), genweight)
                        h['hMuTau_highMt_dRj'].Fill(j.DeltaR(mu+tau), genweight)
                        h['hMuTau_highMt_TauPt'].Fill(tau.Pt(), genweight)
                        h['hMuTau_highMt_MuonPt'].Fill(mu.Pt(), genweight)
                        h['hMuTau_highMt_MetPt'].Fill(met_pt, genweight)
                        h['hMuTau_highMt_JetPt'].Fill(j.Pt(), genweight)
                        h['hMuTau_highMt'].Fill((mu+tau).M(), genweight)
                        h['hMuTau_highMt_Nj'].Fill(len(s_j), genweight)
                        h['hMuTau_highMt_Nbj'].Fill(len(s_b), genweight)

                        if met_pt > 100 : #highMt, highMET
                            h['hMuTau_highMt_highMET_dRl'].Fill(mu.DeltaR(tau), genweight)
                            h['hMuTau_highMt_highMET_dRj'].Fill(j.DeltaR(mu+tau), genweight)
                            h['hMuTau_highMt_highMET_TauPt'].Fill(tau.Pt(), genweight)
                            h['hMuTau_highMt_highMET_MuonPt'].Fill(mu.Pt(), genweight)
                            h['hMuTau_highMt_highMET_MetPt'].Fill(met_pt, genweight)
                            h['hMuTau_highMt_highMET_JetPt'].Fill(j.Pt(), genweight)
                            h['hMuTau_highMt_highMET'].Fill((mu+tau).M(), genweight)
                            h['hMuTau_highMt_highMET_Nj'].Fill(len(s_j), genweight)

                        if mu.DeltaR(tau) < 0.4 and mu.DeltaR(j) > 0.8 and tau.DeltaR(j) > 0.8 : #dR cut                               
                            h['hMuTau_highMt_dRcut_TauPt'].Fill(tau.Pt(), genweight)
                            h['hMuTau_highMt_dRcut_MuonPt'].Fill(mu.Pt(), genweight)
                            h['hMuTau_highMt_dRcut_MetPt'].Fill(met_pt, genweight)
                            h['hMuTau_highMt_dRcut_JetPt'].Fill(j.Pt(), genweight)
                            h['hMuTau_highMt_dRcut'].Fill((mu+tau).M(), genweight)
                            h['hMuTau_highMt_dRcut_Nj'].Fill(len(s_j), genweight)
                            h['hMuTau_highMt_dRcut_Nbj'].Fill(len(s_b), genweight)

                            if met_pt > 100 : #MET cut
                                h['hMuTau_highMt_dRcut_highMET_MuonPt'].Fill(mu.Pt(), genweight)
                                h['hMuTau_highMt_dRcut_highMET_TauPt'].Fill(tau.Pt(), genweight)
                                h['hMuTau_highMt_dRcut_highMET_JetPt'].Fill(j.Pt(), genweight)
                                h['hMuTau_highMt_dRcut_highMET'].Fill((mu+tau).M(), genweight)
                                h['hMuTau_highMt_dRcut_highMET_Nj'].Fill(len(s_j), genweight)
                                h['hMuTau_highMt_dRcut_highMET_Nbj'].Fill(len(s_b), genweight)

                                if mtau[0].decaymode == 0 : h['hMuTau_highMt_dRcut_highMET_TauPt0'].Fill(tau.Pt(), genweight)
                                if mtau[0].decaymode == 1 : h['hMuTau_highMt_dRcut_highMET_TauPt1'].Fill(tau.Pt(), genweight)
                                if mtau[0].decaymode == 10 : h['hMuTau_highMt_dRcut_highMET_TauPt10'].Fill(tau.Pt(), genweight)



                if lmuon[0].charge*mtau[0].charge > 0 :
                    h['hMuTau_SS_Nbj'].Fill(len(s_b), genweight)

                if lmuon[0].charge*mtau[0].charge > 0 and len(s_b) == 0: #SS
                    h['hMuTau_SS_TauPt'].Fill(tau.Pt(), genweight)
                    h['hMuTau_SS_MuonPt'].Fill(mu.Pt(), genweight)
                    h['hMuTau_SS_MetPt'].Fill(met_pt, genweight)
                    h['hMuTau_SS_JetPt'].Fill(j.Pt(), genweight)
                    h['hMuTau_SS'].Fill((mu+tau).M(), genweight)
                    h['hMuTau_SS_Nj'].Fill(len(s_j), genweight)
                    h['hMuTau_SS_Nbj'].Fill(len(s_b), genweight)
                    h['hMuTau_SS_Mt'].Fill(Mt(mu,m), genweight)
                    h['hMuTau_SS_dRl'].Fill(mu.DeltaR(tau), genweight)
                    h['hMuTau_SS_dRj'].Fill(j.DeltaR(mu+tau), genweight)

                    if met_pt < 100 :
                        h['hMuTau_SS_lowMET_TauPt'].Fill(tau.Pt(), genweight)
                        h['hMuTau_SS_lowMET_MuonPt'].Fill(mu.Pt(), genweight)
                        h['hMuTau_SS_lowMET_MetPt'].Fill(met_pt, genweight)
                        h['hMuTau_SS_lowMET_JetPt'].Fill(j.Pt(), genweight)
                        h['hMuTau_SS_lowMET'].Fill((mu+tau).M(), genweight)
                        h['hMuTau_SS_lowMET_Nj'].Fill(len(s_j), genweight)
                        h['hMuTau_SS_lowMET_Nbj'].Fill(len(s_b), genweight)
                        h['hMuTau_SS_lowMET_Mt'].Fill(Mt(mu,m), genweight)
                        h['hMuTau_SS_lowMET_dRl'].Fill(mu.DeltaR(tau), genweight)
                        h['hMuTau_SS_lowMET_dRj'].Fill(j.DeltaR(mu+tau), genweight)

                        if Mt(mu,m) > 50 : #highMt          
                            h['hMuTau_SS_lowMET_highMt_TauPt'].Fill(tau.Pt(), genweight)
                            h['hMuTau_SS_lowMET_highMt_MuonPt'].Fill(mu.Pt(), genweight)
                            h['hMuTau_SS_lowMET_highMt_MetPt'].Fill(met_pt, genweight)
                            h['hMuTau_SS_lowMET_highMt_JetPt'].Fill(j.Pt(), genweight)
                            h['hMuTau_SS_lowMET_highMt'].Fill((mu+tau).M(), genweight)
                            h['hMuTau_SS_lowMET_highMt_Nj'].Fill(len(s_j), genweight)
                            h['hMuTau_SS_lowMET_highMt_Nbj'].Fill(len(s_b), genweight)
                            h['hMuTau_SS_lowMET_highMt_dRl'].Fill(mu.DeltaR(tau), genweight)
                            h['hMuTau_SS_lowMET_highMt_dRj'].Fill(j.DeltaR(mu+tau), genweight)

                        if Mt(mu,m) < 50 : #lowMt   
                            h['hMuTau_SS_lowMET_lowMt_TauPt'].Fill(tau.Pt(), genweight)
                            h['hMuTau_SS_lowMET_lowMt_MuonPt'].Fill(mu.Pt(), genweight)
                            h['hMuTau_SS_lowMET_lowMt_MetPt'].Fill(met_pt, genweight)
                            h['hMuTau_SS_lowMET_lowMt_JetPt'].Fill(j.Pt(), genweight)
                            h['hMuTau_SS_lowMET_lowMt'].Fill((mu+tau).M(), genweight)
                            h['hMuTau_SS_lowMET_lowMt_Nj'].Fill(len(s_j), genweight)
                            h['hMuTau_SS_lowMET_lowMt_Nbj'].Fill(len(s_b), genweight)
                            h['hMuTau_SS_lowMET_lowMt_dRl'].Fill(mu.DeltaR(tau), genweight)
                            h['hMuTau_SS_lowMET_lowMt_dRj'].Fill(j.DeltaR(mu+tau), genweight)

                        if mu.DeltaR(tau) < 0.4 and mu.DeltaR(j) > 0.8 and tau.DeltaR(j) > 0.8 :
                            h['hMuTau_SS_lowMET_dRcut_TauPt'].Fill(tau.Pt(), genweight)
                            h['hMuTau_SS_lowMET_dRcut_MuonPt'].Fill(mu.Pt(), genweight)
                            h['hMuTau_SS_lowMET_dRcut_MetPt'].Fill(met_pt, genweight)
                            h['hMuTau_SS_lowMET_dRcut_JetPt'].Fill(j.Pt(), genweight)
                            h['hMuTau_SS_lowMET_dRcut'].Fill((mu+tau).M(), genweight)
                            h['hMuTau_SS_lowMET_dRcut_Nj'].Fill(len(s_j), genweight)
                            h['hMuTau_SS_lowMET_dRcut_Nbj'].Fill(len(s_b), genweight)
                            h['hMuTau_SS_lowMET_dRcut_Mt'].Fill(Mt(mu,m), genweight)

                            if Mt(mu,m) > 50 : #highMt
                                h['hMuTau_SS_lowMET_dRcut_highMt_TauPt'].Fill(tau.Pt(), genweight)
                                h['hMuTau_SS_lowMET_dRcut_highMt_MuonPt'].Fill(mu.Pt(), genweight)
                                h['hMuTau_SS_lowMET_dRcut_highMt_MetPt'].Fill(met_pt, genweight)
                                h['hMuTau_SS_lowMET_dRcut_highMt_JetPt'].Fill(j.Pt(), genweight)
                                h['hMuTau_SS_lowMET_dRcut_highMt'].Fill((mu+tau).M(), genweight)
                                h['hMuTau_SS_lowMET_dRcut_highMt_Nj'].Fill(len(s_j), genweight)
                                h['hMuTau_SS_lowMET_dRcut_highMt_Nbj'].Fill(len(s_b), genweight)

                                if mtau[0].decaymode == 0 : h['hMuTau_SS_lowMET_dRcut_highMt_TauPt0'].Fill(tau.Pt(), genweight)
                                if mtau[0].decaymode == 1 : h['hMuTau_SS_lowMET_dRcut_highMt_TauPt1'].Fill(tau.Pt(), genweight)
                                if mtau[0].decaymode == 10 : h['hMuTau_SS_lowMET_dRcut_highMt_TauPt10'].Fill(tau.Pt(), genweight)

                                h['hMuTau_SS_lowMET_dRcut_highMt_TauPtJetPt'].Fill(tau.Pt(), j.Pt(), genweight)
                                h['hMuTau_SS_lowMET_dRcut_highMt_TauPtMuonPt'].Fill(tau.Pt(), mu.Pt(), genweight)
                                h['hMuTau_SS_lowMET_dRcut_highMt_TauPtMetPt'].Fill(tau.Pt(), met_pt, genweight)

                                if len(s_mu) > 1:
                                    h['hMuTau_SS_lowMET_dRcut_highMt_TauPtMuMuMass'].Fill(tau.Pt(), (mu+mu2).M(), genweight)

                                if len(s_j) > 1:
                                    h['hMuTau_SS_lowMET_dRcut_highMt_TauPtdRj'].Fill(tau.Pt(), tau.DeltaR(j), genweight)
                                    h['hMuTau_SS_lowMET_dRcut_highMt_TauPtdRj2'].Fill(tau.Pt(), tau.DeltaR(j2), genweight)

                            if Mt(mu,m) < 50 : #lowMt
                                h['hMuTau_SS_lowMET_dRcut_lowMt_TauPt'].Fill(tau.Pt(), genweight)
                                h['hMuTau_SS_lowMET_dRcut_lowMt_MuonPt'].Fill(mu.Pt(), genweight)
                                h['hMuTau_SS_lowMET_dRcut_lowMt_MetPt'].Fill(met_pt, genweight)
                                h['hMuTau_SS_lowMET_dRcut_lowMt_JetPt'].Fill(j.Pt(), genweight)
                                h['hMuTau_SS_lowMET_dRcut_lowMt'].Fill((mu+tau).M(), genweight)
                                h['hMuTau_SS_lowMET_dRcut_lowMt_Nj'].Fill(len(s_j), genweight)
                                h['hMuTau_SS_lowMET_dRcut_lowMt_Nbj'].Fill(len(s_b), genweight)

                                if mtau[0].decaymode == 0 : h['hMuTau_SS_lowMET_dRcut_lowMt_TauPt0'].Fill(tau.Pt(), genweight)
                                if mtau[0].decaymode == 1 : h['hMuTau_SS_lowMET_dRcut_lowMt_TauPt1'].Fill(tau.Pt(), genweight)
                                if mtau[0].decaymode == 10 : h['hMuTau_SS_lowMET_dRcut_lowMt_TauPt10'].Fill(tau.Pt(), genweight)

                                h['hMuTau_SS_lowMET_dRcut_lowMt_TauPtJetPt'].Fill(tau.Pt(), j.Pt(), genweight)
                                h['hMuTau_SS_lowMET_dRcut_lowMt_TauPtMuonPt'].Fill(tau.Pt(), mu.Pt(), genweight)
                                h['hMuTau_SS_lowMET_dRcut_lowMt_TauPtMetPt'].Fill(tau.Pt(), met_pt, genweight)

                                if len(s_mu) > 1:
                                    h['hMuTau_SS_lowMET_dRcut_lowMt_TauPtMuMuMass'].Fill(tau.Pt(), (mu+mu2).M(), genweight)

                                if len(s_j) > 1:
                                    h['hMuTau_SS_lowMET_dRcut_lowMt_TauPtdRj'].Fill(tau.Pt(), tau.DeltaR(j), genweight)
                                    h['hMuTau_SS_lowMET_dRcut_lowMt_TauPtdRj2'].Fill(tau.Pt(), tau.DeltaR(j2), genweight)

                    if mu.DeltaR(tau) < 0.4 and mu.DeltaR(j) > 0.8 and tau.DeltaR(j) > 0.8 : #dR cut
                        h['hMuTau_SS_dRcut_TauPt'].Fill(tau.Pt(), genweight)
                        h['hMuTau_SS_dRcut_MuonPt'].Fill(mu.Pt(), genweight)
                        h['hMuTau_SS_dRcut_MetPt'].Fill(met_pt, genweight)
                        h['hMuTau_SS_dRcut_JetPt'].Fill(j.Pt(), genweight)
                        h['hMuTau_SS_dRcut'].Fill((mu+tau).M(), genweight)
                        h['hMuTau_SS_dRcut_Nj'].Fill(len(s_j), genweight)
                        h['hMuTau_SS_dRcut_Nbj'].Fill(len(s_b), genweight)
                        h['hMuTau_SS_dRcut_Mt'].Fill(Mt(mu,m), genweight)

                        if met_pt > 100 : #Met cut
                            h['hMuTau_SS_dRcut_highMET_MuonPt'].Fill(mu.Pt(), genweight)
                            h['hMuTau_SS_dRcut_highMET_MetPt'].Fill(met_pt, genweight)
                            h['hMuTau_SS_dRcut_highMET_TauPt'].Fill(tau.Pt(), genweight)
                            h['hMuTau_SS_dRcut_highMET_JetPt'].Fill(j.Pt(), genweight)
                            h['hMuTau_SS_dRcut_highMET'].Fill((mu+tau).M(), genweight)
                            h['hMuTau_SS_dRcut_highMET_Nj'].Fill(len(s_j), genweight)
                            h['hMuTau_SS_dRcut_highMET_Nbj'].Fill(len(s_b), genweight)
                            h['hMuTau_SS_dRcut_highMET_Mt'].Fill(Mt(mu,m), genweight)

                            if Mt(mu,m) > 50 : #high mT
                                h['hMuTau_SS_dRcut_highMET_highMt_TauPt'].Fill(tau.Pt(), genweight)
                                h['hMuTau_SS_dRcut_highMET_highMt_MuonPt'].Fill(mu.Pt(), genweight)
                                h['hMuTau_SS_dRcut_highMET_highMt_MetPt'].Fill(met_pt, genweight)
                                h['hMuTau_SS_dRcut_highMET_highMt_JetPt'].Fill(j.Pt(), genweight)
                                h['hMuTau_SS_dRcut_highMET_highMt'].Fill((mu+tau).M(), genweight)
                                h['hMuTau_SS_dRcut_highMET_highMt_Nj'].Fill(len(s_j), genweight)
                                h['hMuTau_SS_dRcut_highMET_highMt_Nbj'].Fill(len(s_b), genweight)

                                if mtau[0].decaymode == 0 : h['hMuTau_SS_dRcut_highMET_highMt_TauPt0'].Fill(tau.Pt(), genweight)
                                if mtau[0].decaymode == 1 : h['hMuTau_SS_dRcut_highMET_highMt_TauPt1'].Fill(tau.Pt(), genweight)
                                if mtau[0].decaymode == 10 : h['hMuTau_SS_dRcut_highMET_highMt_TauPt10'].Fill(tau.Pt(), genweight)

                                h['hMuTau_SS_dRcut_highMET_highMt_TauPtJetPt'].Fill(tau.Pt(), j.Pt(), genweight)
                                h['hMuTau_SS_dRcut_highMET_highMt_TauPtMuonPt'].Fill(tau.Pt(), mu.Pt(), genweight)
                                h['hMuTau_SS_dRcut_highMET_highMt_TauPtMetPt'].Fill(tau.Pt(), met_pt, genweight)

                                if len(s_mu) > 1:
                                    h['hMuTau_SS_dRcut_highMET_highMt_TauPtMuMuMass'].Fill(tau.Pt(), (mu+mu2).M(), genweight)

                                if len(s_j) > 1:
                                    h['hMuTau_SS_dRcut_highMET_highMt_TauPtdRj'].Fill(tau.Pt(), tau.DeltaR(j), genweight)
                                    h['hMuTau_SS_dRcut_highMET_highMt_TauPtdRj2'].Fill(tau.Pt(), tau.DeltaR(j2), genweight)


                            if Mt(mu,m) < 50 : #mT cut #DR3
                                h['hMuTau_SS_dRcut_highMET_lowMt_TauPt'].Fill(tau.Pt(), genweight)
                                h['hMuTau_SS_dRcut_highMET_lowMt_MuonPt'].Fill(mu.Pt(), genweight)
                                h['hMuTau_SS_dRcut_highMET_lowMt_MetPt'].Fill(met_pt, genweight)
                                h['hMuTau_SS_dRcut_highMET_lowMt_JetPt'].Fill(j.Pt(), genweight)
                                h['hMuTau_SS_dRcut_highMET_lowMt'].Fill((mu+tau).M(), genweight)
                                h['hMuTau_SS_dRcut_highMET_lowMt_Nj'].Fill(len(s_j), genweight)
                                h['hMuTau_SS_dRcut_highMET_lowMt_Nbj'].Fill(len(s_b), genweight)

                                if mtau[0].decaymode == 0 : h['hMuTau_SS_dRcut_highMET_lowMt_TauPt0'].Fill(tau.Pt(), genweight)
                                if mtau[0].decaymode == 1 : h['hMuTau_SS_dRcut_highMET_lowMt_TauPt1'].Fill(tau.Pt(), genweight)
                                if mtau[0].decaymode == 10 : h['hMuTau_SS_dRcut_highMET_lowMt_TauPt10'].Fill(tau.Pt(), genweight)

                                h['hMuTau_SS_dRcut_highMET_lowMt_TauPtJetPt'].Fill(tau.Pt(), j.Pt(), genweight)
                                h['hMuTau_SS_dRcut_highMET_lowMt_TauPtMuonPt'].Fill(tau.Pt(), mu.Pt(), genweight)
                                h['hMuTau_SS_dRcut_highMET_lowMt_TauPtMetPt'].Fill(tau.Pt(), met_pt, genweight)
                                
                                if len(s_mu) > 1:
                                    h['hMuTau_SS_dRcut_highMET_lowMt_TauPtMuMuMass'].Fill(tau.Pt(), (mu+mu2).M(), genweight)

                                if len(s_j) > 1:
                                    h['hMuTau_SS_dRcut_highMET_lowMt_TauPtdRj'].Fill(tau.Pt(), tau.DeltaR(j), genweight)
                                    h['hMuTau_SS_dRcut_highMET_lowMt_TauPtdRj2'].Fill(tau.Pt(), tau.DeltaR(j2), genweight)


       if isData == 0 or isAltered == 1:
            if isJetHTEvent == 1 or isSingleMuonEvent == 1 :
                if j.Pt() > 100 and (mu+tau).M() > 1:
                    if lmuon[0].charge*mtau[0].charge < 0 and len(s_b) == 0 :
                        h['hMuTau_SR_Events'].Fill(1, genweight)
                        h['hMuTau_SR_dRl'].Fill(mu.DeltaR(tau), genweight)
                        h['hMuTau_SR_dRj'].Fill(j.DeltaR(mu+tau), genweight)
                        h['hMuTau_SR_TauPt'].Fill(tau.Pt(), genweight)
                        h['hMuTau_SR_MuonPt'].Fill(mu.Pt(), genweight)
                        h['hMuTau_SR_MetPt'].Fill(met_pt, genweight)
                        h['hMuTau_SR_JetPt'].Fill(j.Pt(), genweight)
                        h['hMuTau_SR'].Fill((mu+tau).M(), genweight)
                        h['hMuTau_SR_Nj'].Fill(len(s_j), genweight)
                        h['hMuTau_SR_Nbj'].Fill(len(s_b), genweight)
                        h['hMuTau_SR_Mt'].Fill(Mt(mu,m), genweight)

                        if mu.DeltaR(tau) < 0.4 and mu.DeltaR(j) > 0.8 and tau.DeltaR(j) > 0.8 :
                            h['hMuTau_SR_Events'].Fill(2, genweight)
                            h['hMuTau_SR_dRcut_dRl'].Fill(mu.DeltaR(tau), genweight)
                            h['hMuTau_SR_dRcut_dRj'].Fill(j.DeltaR(mu+tau), genweight)
                            h['hMuTau_SR_dRcut_TauPt'].Fill(tau.Pt(), genweight)
                            h['hMuTau_SR_dRcut_MuonPt'].Fill(mu.Pt(), genweight)
                            h['hMuTau_SR_dRcut_MetPt'].Fill(met_pt, genweight)
                            h['hMuTau_SR_dRcut_JetPt'].Fill(j.Pt(), genweight)
                            h['hMuTau_SR_dRcut'].Fill((mu+tau).M(), genweight)
                            h['hMuTau_SR_dRcut_Nj'].Fill(len(s_j), genweight)
                            h['hMuTau_SR_dRcut_Nbj'].Fill(len(s_b), genweight)
                            h['hMuTau_SR_dRcut_Mt'].Fill(Mt(mu,m), genweight)

                            if met_pt > 100 :
                                h['hMuTau_SR_Events'].Fill(3, genweight)
                                h['hMuTau_SR_dRcut_highMET_dRl'].Fill(mu.DeltaR(tau), genweight)
                                h['hMuTau_SR_dRcut_highMET_dRj'].Fill(j.DeltaR(mu+tau), genweight)
                                h['hMuTau_SR_dRcut_highMET_TauPt'].Fill(tau.Pt(), genweight)
                                h['hMuTau_SR_dRcut_highMET_MuonPt'].Fill(mu.Pt(), genweight)
                                h['hMuTau_SR_dRcut_highMET_MetPt'].Fill(met_pt, genweight)
                                h['hMuTau_SR_dRcut_highMET_JetPt'].Fill(j.Pt(), genweight)
                                h['hMuTau_SR_dRcut_highMET'].Fill((mu+tau).M(), genweight)
                                h['hMuTau_SR_dRcut_highMET_Nj'].Fill(len(s_j), genweight)
                                h['hMuTau_SR_dRcut_highMET_Nbj'].Fill(len(s_b), genweight)
                                h['hMuTau_SR_dRcut_highMET_Mt'].Fill(Mt(mu,m), genweight)

                                if Mt(mu,m) < 50 :
                                    h['hMuTau_SR_Events'].Fill(4, genweight)
                                    h['hMuTau_SR_dRcut_highMET_lowMt_dRl'].Fill(mu.DeltaR(tau), genweight)
                                    h['hMuTau_SR_dRcut_highMET_lowMt_dRj'].Fill(j.DeltaR(mu+tau), genweight)
                                    h['hMuTau_SR_dRcut_highMET_lowMt_TauPt'].Fill(tau.Pt(), genweight)
                                    h['hMuTau_SR_dRcut_highMET_lowMt_MuonPt'].Fill(mu.Pt(), genweight)
                                    h['hMuTau_SR_dRcut_highMET_lowMt_MetPt'].Fill(met_pt, genweight)
                                    h['hMuTau_SR_dRcut_highMET_lowMt_JetPt'].Fill(j.Pt(), genweight)
                                    h['hMuTau_SR_dRcut_highMET_lowMt'].Fill((mu+tau).M(), genweight)
                                    h['hMuTau_SR_dRcut_highMET_lowMt_Nj'].Fill(len(s_j), genweight)
                                    h['hMuTau_SR_dRcut_highMET_lowMt_Nbj'].Fill(len(s_b), genweight)
                                    h['hMuTau_SR_dRcut_highMET_lowMt_Mt'].Fill(Mt(mu,m), genweight)

                                    if mtau[0].decaymode == 0 : h['hMuTau_SR_dRcut_highMET_lowMt_TauPt0'].Fill(tau.Pt(), genweight)
                                    if mtau[0].decaymode == 1 : h['hMuTau_SR_dRcut_highMET_lowMt_TauPt1'].Fill(tau.Pt(), genweight)
                                    if mtau[0].decaymode == 10 : h['hMuTau_SR_dRcut_highMET_lowMt_TauPt10'].Fill(tau.Pt(), genweight)

                                    if isJetHTEvent == 1 : h['hMuTau_SR_Trigger'].Fill(7, genweight)
                                    h['hMuTau_SR_Inclusive_Trigger'].Fill(7, genweight)

            if j.Pt() > 100 and (mu+tau).M() > 1:
                if lmuon[0].charge*mtau[0].charge < 0 and len(s_b) == 0 :
                    if mu.DeltaR(tau) < 0.4 and mu.DeltaR(j) > 0.8 and tau.DeltaR(j) > 0.8 :
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



       return isMuTau



inputFileNames=open(inputFileList, 'r')
for inputFileName in inputFileNames:
    inputFileName=inputFileName.replace("\n","")
    print(inputFileName.replace("\n",""))

    fchain.Add(inputFileName)
    chain2.Add(inputFileName)
    if isData == 0:
        chain3.Add(inputFileName)


fchain.AddFriend(chain2)
if isData == 0:
    fchain.AddFriend(chain3)

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

   h['hEvents'].Fill(0.5, 1)

   h['hMetPt'].Fill(met_pt)

   genweight = 1

   if isData == 0:
       weight = fchain.GetBranch("lumiInfo")
       genweight = weight.GetLeaf('weight').GetValue()

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

   s_j = []
   s_b = []
   s_e = []
   s_mu = []
   s_isomu = []
   s_isoe = []

   unclean = []
   eclean = []
   mclean = []

   boosted = []
   mclean_altered = []

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
           if tau.mvaid >= 1 and tau.mvaid < 4 :
               mclean_altered+=[tau]
               h['hMuCleanedPt_altered'].Fill(tau.pt, genweight)
           if tau.mvaid >= 4 :
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
             h['hJetPt'].Fill(jet.pt, genweight)
             h['hBtag'].Fill(jet.deepjet)
             s_j+=[jet]
             if jet.deepjet > 0.7476:
                 h['hBJetPt'].Fill(jet.pt, genweight)
                 s_b+=[jet]

   if muons.size()>0:
      for i in range(muons.size()):
         muon = muons.at(i)
         if muon.id >= 1:
             h['hMuonPt'].Fill(muon.pt, genweight) 
             s_mu+=[muon]
             if muon.iso < 0.25:
                 h['hIsoMuonPt'].Fill(muon.pt, genweight)
                 s_isomu+=[muon]

   if electrons.size()>0:
      for i in range(electrons.size()):
         electron = electrons.at(i)
         if electron.id >= 1 :
             h['hElectronPt'].Fill(electron.pt, genweight)
             s_e+=[electron]
             if electron.iso >= 1:
                 h['hIsoElectronPt'].Fill(electron.pt, genweight)
                 s_isoe+=[electron]

   s_j.sort(key=lambda x: x.pt, reverse=True)
   s_e.sort(key=lambda x: x.pt, reverse=True)
   s_isoe.sort(key=lambda x: x.pt, reverse=True)
   s_mu.sort(key=lambda x: x.pt, reverse=True)
   s_isomu.sort(key=lambda x: x.pt, reverse=True)
   unclean.sort(key=lambda x: x.pt, reverse=True)
   eclean.sort(key=lambda x: x.pt, reverse=True)
   mclean.sort(key=lambda x: x.pt, reverse=True)
   mclean_altered.sort(key=lambda x: x.pt, reverse=True)

   boosted.sort(key=lambda x: x.pt, reverse=True)


   if len(s_isomu) > 1 and len(s_j) > 0 and s_isomu[0].charge*s_isomu[1].charge < 0 : 
       if MuMu_Channel(s_isomu) == 1: continue

   if len(s_isomu) > 0 and len(s_isoe) > 0 and len(s_j) > 0 and len(s_b) == 0 and s_isoe[0].charge*s_isomu[0].charge < 0 : 
       if EMu_Channel(s_isoe,s_isomu) == 1: continue

   if isAltered == 1 :
       if len(s_mu) > 0 and len(mclean_altered) > 0 and len(s_j) > 0 : 
           MuTau_Channel(mclean_altered, s_mu)

   if isAltered == 0 :
       if len(s_mu) > 0 and len(mclean) > 0 and len(s_j) > 0 : 
           MuTau_Channel(mclean, s_mu)

#   if len(s_isoe) > 1 and len(s_j) > 0 and len(s_b) == 0 and s_isoe[0].charge*s_isoe[1].charge < 0 : 
#       if EE_Channel(s_isoe) == 1: continue



out.cd()

for key in h.keys():
    h[key].Write()

out.Close()

print("--- %s seconds ---" % (time.time() - start_time))
