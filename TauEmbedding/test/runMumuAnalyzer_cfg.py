import os
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

isMC = True
whichTrigger="DoubleMu"
#whichTrigger="SingleMu"
#whichTrigger="JetHT"

prefix = 'root://xrootd.unl.edu/'

if not isMC:
    if whichTrigger == "DoubleMu":
        inputFiles = '/store/data/Run2016B/DoubleMuon/MINIAOD/17Jul2018_ver1-v1/50000/3AA00BFC-258B-E811-97ED-0090FAA58194.root'
    elif whichTrigger == "SingleMu":
        inputFiles = '/store/data/Run2016B/SingleMuon/MINIAOD/17Jul2018_ver1-v1/80000/306DAB6C-068C-E811-9E30-0242AC1C0501.root'
    elif whichTrigger == "JetHT":
        inputFiles = '/store/data/Run2016B/JetHT/MINIAOD/17Jul2018_ver2-v2/100000/00323A27-15B8-E811-88B5-E0071B6C9DB0.root'
else:
    inputFiles = '/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-5to50_HT-70to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/100000/026F8D54-5AEA-E811-91CC-24BE05CE2D41.root'

if not isMC:
    outputFile = "outMuonAnalyzer"+whichTrigger+".root"
else:
    outputFile = "outMuonAnalyzerMC"+whichTrigger+".root"
maxEvents = -1

process = cms.Process('Analysis',eras.Run2_2017,eras.run2_miniAOD_94XFall17)

# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(maxEvents)
)

# Input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(prefix+inputFiles),
    secondaryFileNames = cms.untracked.vstring(),
)

# Output
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outputFile),
                                   closeFileFast = cms.untracked.bool(True)
)

process.load('BoostedDiTau.TauEmbedding.SelectingProcedure_cff')
if not isMC:
    if whichTrigger == "DoubleMu":
        process.p = cms.Path(process.studyPatMuonsZmumuSelectionDoubleMu)
    elif whichTrigger == "SingleMu":
        process.p = cms.Path(process.studyPatMuonsZmumuSelectionSingleMu)
        process.ZmumuCandidates.cut = cms.string('charge = 0 & max(daughter(0).pt, daughter(1).pt) > 24 & mass > 1 & daughter(0).isGlobalMuon & daughter(1).isGlobalMuon')
        process.ZmumuCandidates.decay=cms.string("patMuonsAfterLooseIDSingleMu@+ patMuonsAfterLooseIDSingleMu@-")
    elif whichTrigger == "JetHT":
        process.p = cms.Path(process.studyPatMuonsZmumuSelectionJetHT)
        process.ZmumuCandidates.cut = cms.string('charge = 0 & mass > 1 & daughter(0).isGlobalMuon & daughter(1).isGlobalMuon')
        process.ZmumuCandidates.decay=cms.string("patMuonsAfterLooseIDJetHT@+ patMuonsAfterLooseIDJetHT@-")
else:
    process.p = cms.Path(process.studyPatMuonsZmumuSelectionMC)
    
process.schedule = cms.Schedule(process.p)
