import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.PAT_cff import *

from PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi import patMuons
from HLTrigger.HLTfilters.triggerResultsFilter_cfi import *
from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector

## Trigger requirements
doubleMuonHLTTrigger = cms.EDFilter("TriggerResultsFilter",
    hltResults = cms.InputTag("TriggerResults","","HLT"),
    l1tResults = cms.InputTag(""),
    throw = cms.bool(False),
    triggerConditions = cms.vstring("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v* OR HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*")
)

singleMuonHLTTrigger = cms.EDFilter("TriggerResultsFilter",
    hltResults = cms.InputTag("TriggerResults","","HLT"),
    l1tResults = cms.InputTag(""),
    throw = cms.bool(False),
    triggerConditions = cms.vstring("HLT_IsoMu24_v* OR HLT_IsoTkMu24_v*")
)

jetHTHLTTrigger = cms.EDFilter("TriggerResultsFilter",
    hltResults = cms.InputTag("TriggerResults","","HLT"),
    l1tResults = cms.InputTag(""),
    throw = cms.bool(False),
    triggerConditions = cms.vstring("HLT_PFJet450_v* OR HLT_PFHT900_v*")
)

## hltTriggerMC = cms.EDFilter("TriggerResultsFilter",
##     hltResults = cms.InputTag("TriggerResults","","HLT"),
##     l1tResults = cms.InputTag(""),
##     throw = cms.bool(False),
##     triggerConditions = cms.vstring("HLT_PFJet450_v* OR HLT_PFHT900_v* OR HLT_IsoMu24_v* OR HLT_IsoTkMu24_v* OR HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v* OR HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*")
## )

## ## Jet selection
## patJetsAfterKinCuts = cms.EDFilter('PATJetSelector',
##                                    src = cms.InputTag('slimmedJets'),
##                                    cut = cms.string('pt > 500.0 && abs(eta)<2.5'),
##                                    filter = cms.bool(True)
## )
## 
## patJetsAfterLooseID = cms.EDFilter("PFJetIDSelectionFunctorFilter",
##                                    filterParams = cms.PSet(
##                                        version = cms.string('WINTER16'),
##                                        quality = cms.string('LOOSE')
##                                    ),
##                                    src = cms.InputTag("patJetsAfterKinCuts"),
##                                    filter = cms.bool(True)
## )

## Muon selection
patMuonsAfterKinCutsDoubleMu = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("slimmedMuons"),
    cut = cms.string("pt > 8 && abs(eta) < 2.5"),
    filter = cms.bool(True)
)

patMuonsAfterKinCutsSingleMu = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("slimmedMuons"),
    cut = cms.string("pt > 3 && abs(eta) < 2.5"),
    filter = cms.bool(True)
)

patMuonsAfterKinCutsJetHT = patMuonsAfterKinCutsSingleMu.clone()

## patMuonsAfterKinCutsMC = patMuonsAfterKinCutsSingleMu.clone()

# For impact parameter (w.r.t. to PV) requirements, a vector collection is needed, therefore only dB < 0.2 required.
# The default requirements (in C++):
# 1) fabs(recoMu.muonBestTrack()->dxy(vertex->position())) < 0.2 ----> similar to dB < 0.2
# 2) fabs(recoMu.muonBestTrack()->dz(vertex->position())) < 0.5
## patMuonsAfterTightID = cms.EDFilter("PATMuonSelector",
##     src = cms.InputTag("patMuonsAfterKinCuts"),
##     cut = cms.string(
##     "isPFMuon && isGlobalMuon"
##     " && muonID('GlobalMuonPromptTight')"
##     " && numberOfMatchedStations > 1"
##     " && innerTrack.hitPattern.trackerLayersWithMeasurement > 5"
##     " && innerTrack.hitPattern.numberOfValidPixelHits > 0"
##     " && dB < 0.2"
##     ),
##     filter = cms.bool(True)
## )
## 
## patMuonsAfterMediumID = cms.EDFilter("PATMuonSelector",
##     src = cms.InputTag("patMuonsAfterKinCuts"),
##     cut = cms.string("isMediumMuon"),
##     filter = cms.bool(True)
## )

patMuonsAfterLooseIDDoubleMu = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsAfterKinCutsDoubleMu"),
    cut = cms.string("isLooseMuon"),
    filter = cms.bool(True)
)

patMuonsAfterLooseIDLooseIsoDoubleMu = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsAfterKinCutsDoubleMu"),
    cut = cms.string("isLooseMuon && (pfIsolationR04().sumChargedHadronPt + max(0., pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - 0.5*pfIsolationR04().sumPUPt))/pt < 0.25"),
    filter = cms.bool(True)
)

patMuonsAfterLooseIDLooseIsoDBDoubleMu = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsAfterKinCutsDoubleMu"),
    cut = cms.string("isLooseMuon && (pfIsolationR04().sumChargedHadronPt + max(0., pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - 0.5*pfIsolationR04().sumPUPt))/pt < 0.25 && dB < 0.2"),
    filter = cms.bool(True)
)

patMuonsAfterLooseIDSingleMu = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsAfterKinCutsSingleMu"),
    cut = cms.string("isLooseMuon"),
    filter = cms.bool(True)
)

patMuonsAfterLooseIDJetHT = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsAfterKinCutsJetHT"),
    cut = cms.string("isLooseMuon"),
    filter = cms.bool(True)
)

## patMuonsAfterLooseIDMC = cms.EDFilter("PATMuonSelector",
##     src = cms.InputTag("patMuonsAfterKinCutsMC"),
##     cut = cms.string("isLooseMuon"),
##     filter = cms.bool(True)
## )

#patMuonsAfterID = patMuonsAfterLooseIDDoubleMu.clone()

## Zmumu selection
ZmumuCandidates = cms.EDProducer("CandViewShallowCloneCombiner",
    checkCharge = cms.bool(True),
    # require one of the muons with pT > 17 GeV, and an invariant mass > 1 GeV
    cut = cms.string('charge = 0 & max(daughter(0).pt, daughter(1).pt) > 17 & mass > 70 & daughter(0).isGlobalMuon & daughter(1).isGlobalMuon'),
    decay = cms.string("patMuonsAfterLooseIDDoubleMu@+ patMuonsAfterLooseIDDoubleMu@-")
)

ZmumuCandidatesLooseIso = cms.EDProducer("CandViewShallowCloneCombiner",
    checkCharge = cms.bool(True),
    # require one of the muons with pT > 17 GeV, and an invariant mass > 1 GeV
    cut = cms.string('charge = 0 & max(daughter(0).pt, daughter(1).pt) > 17 & mass > 1 & daughter(0).isGlobalMuon & daughter(1).isGlobalMuon'),
    decay = cms.string("patMuonsAfterLooseIDLooseIsoDoubleMu@+ patMuonsAfterLooseIDLooseIsoDoubleMu@-")
)

ZmumuCandidatesLooseIsoDB = cms.EDProducer("CandViewShallowCloneCombiner",
    checkCharge = cms.bool(True),
    # require one of the muons with pT > 17 GeV, and an invariant mass > 1 GeV
    cut = cms.string('charge = 0 & max(daughter(0).pt, daughter(1).pt) > 17 & mass > 1 & daughter(0).isGlobalMuon & daughter(1).isGlobalMuon'),
    decay = cms.string("patMuonsAfterLooseIDLooseIsoDBDoubleMu@+ patMuonsAfterLooseIDLooseIsoDBDoubleMu@-")
)

ZmumuCandidatesFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("ZmumuCandidates"),
    minNumber = cms.uint32(1)
    #filter = cms.bool(True)
)

ZmumuCandidatesFilterLooseIso = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("ZmumuCandidatesLooseIso"),
    minNumber = cms.uint32(1)
    #filter = cms.bool(True)
)

ZmumuCandidatesFilterLooseIsoDB = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("ZmumuCandidatesLooseIsoDB"),
    minNumber = cms.uint32(1),
    #filter = cms.bool(True)
)

## ZmumuCandidatesMCDoubleMu = cms.EDProducer("CandViewShallowCloneCombiner",
##     checkCharge = cms.bool(True),
##     # require one of the muons with pT > 17 GeV, and an invariant mass > 1 GeV
##     cut = cms.string('charge = 0 & max(daughter(0).pt, daughter(1).pt) > 17 & min(daughter(0).pt, daughter(1).pt) > 8 & mass > 1 & daughter(0).isGlobalMuon & daughter(1).isGlobalMuon'),
##     decay = cms.string("patMuonsAfterLooseIDMC@+ patMuonsAfterLooseIDMC@-")
## )
## 
## ZmumuCandidatesFilterMCDoubleMu = cms.EDFilter("CandViewCountFilter",
##     src = cms.InputTag("ZmumuCandidatesMCDoubleMu"),
##     minNumber = cms.uint32(1),
##     filter = cms.bool(True)
## )
## 
## ZmumuCandidatesMCSingleMu = cms.EDProducer("CandViewShallowCloneCombiner",
##     checkCharge = cms.bool(True),
##     # require one of the muons with pT > 17 GeV, and an invariant mass > 1 GeV
##     cut = cms.string('charge = 0 & max(daughter(0).pt, daughter(1).pt) > 24 & mass > 1 & daughter(0).isGlobalMuon & daughter(1).isGlobalMuon'),
##     decay = cms.string("patMuonsAfterLooseIDMC@+ patMuonsAfterLooseIDMC@-")
## )
## 
## ZmumuCandidatesFilterMCSingleMu = cms.EDFilter("CandViewCountFilter",
##     src = cms.InputTag("ZmumuCandidatesMCSingleMu"),
##     minNumber = cms.uint32(1),
##     filter = cms.bool(True)
## )
## 
## ZmumuCandidatesMCJetHT = cms.EDProducer("CandViewShallowCloneCombiner",
##     checkCharge = cms.bool(True),
##     # require one of the muons with pT > 17 GeV, and an invariant mass > 1 GeV
##     cut = cms.string('charge = 0 & mass > 1 & daughter(0).isGlobalMuon & daughter(1).isGlobalMuon'),
##     decay = cms.string("patMuonsAfterLooseIDMC@+ patMuonsAfterLooseIDMC@-")
## )
## 
## ZmumuCandidatesFilterMCJetHT = cms.EDFilter("CandViewCountFilter",
##     src = cms.InputTag("ZmumuCandidatesMCJetHT"),
##     minNumber = cms.uint32(1),
##     filter = cms.bool(True)
## )


selectedMuonsForEmbedding = cms.EDProducer("MuMuForEmbeddingSelector",
    ZmumuCandidatesCollection = cms.InputTag("ZmumuCandidates")
)
## 
## selectedBoostedMuonsForEmbedding = cms.EDProducer("MuMuForEmbeddingBoostedSelector",
##     ZmumuCandidatesCollection = cms.InputTag("ZmumuCandidates")
## )

## HighPtJetFilter = cms.EDFilter("CandCountFilter",
##                                src = cms.InputTag("patJetsAfterLooseID","Analysis"),
##                                minNumber = cms.uint32(1),
##                                filter = cms.bool(True)
## )

## Analyzer
analyzeMuonsForEmbedding = cms.EDAnalyzer("MuMuForEmbeddingAnalyzer",
                                          ZmumuCandidatesCollection = cms.InputTag("ZmumuCandidates"),
                                          JetCollection = cms.InputTag("slimmedJets"),
                                          METCollection = cms.InputTag("slimmedMETs")
                                          
)

analyzeMuonsForEmbeddingLooseIso = cms.EDAnalyzer("MuMuForEmbeddingAnalyzer",
                                          ZmumuCandidatesCollection = cms.InputTag("ZmumuCandidatesLooseIso"),
                                          JetCollection = cms.InputTag("slimmedJets"),
                                          METCollection = cms.InputTag("slimmedMETs")
                                          
)

analyzeMuonsForEmbeddingLooseIsoDB = cms.EDAnalyzer("MuMuForEmbeddingAnalyzer",
                                          ZmumuCandidatesCollection = cms.InputTag("ZmumuCandidatesLooseIsoDB"),
                                          JetCollection = cms.InputTag("slimmedJets"),
                                          METCollection = cms.InputTag("slimmedMETs")
                                          
)

## analyzeMuonsForEmbeddingMCDoubleMu = cms.EDAnalyzer("MuMuForEmbeddingAnalyzer",
##                                           ZmumuCandidatesCollection = cms.InputTag("ZmumuCandidatesMCDoubleMu"),
##                                           JetCollection = cms.InputTag("slimmedJets"),
##                                           METCollection = cms.InputTag("slimmedMETs")
##                                           
## )
## 
## analyzeMuonsForEmbeddingMCSingleMu = cms.EDAnalyzer("MuMuForEmbeddingAnalyzer",
##                                           ZmumuCandidatesCollection = cms.InputTag("ZmumuCandidatesMCSingleMu"),
##                                           JetCollection = cms.InputTag("slimmedJets"),
##                                           METCollection = cms.InputTag("slimmedMETs")
##                                           
## )
## 
## analyzeMuonsForEmbeddingMCJetHT = cms.EDAnalyzer("MuMuForEmbeddingAnalyzer",
##                                           ZmumuCandidatesCollection = cms.InputTag("ZmumuCandidatesMCJetHT"),
##                                           JetCollection = cms.InputTag("slimmedJets"),
##                                           METCollection = cms.InputTag("slimmedMETs")
##                                           
## )

## analyzeTriggerDoubleMu = cms.EDAnalyzer("MuMuTriggerAnalyzer",
##                                         TriggerResults = cms.InputTag("TriggerResults"),
##                                         ZmumuCandidatesCollection = cms.InputTag("ZmumuCandidatesMCDoubleMu")
## )
## 
## analyzeTriggerSingleMu = cms.EDAnalyzer("MuMuTriggerAnalyzer",
##                                         TriggerResults = cms.InputTag("TriggerResults"),
##                                         ZmumuCandidatesCollection = cms.InputTag("ZmumuCandidatesMCSingleMu")
## )
## 
## analyzeTriggerDoubleMu = cms.EDAnalyzer("MuMuTriggerAnalyzer",
##                                         TriggerResults = cms.InputTag("TriggerResults"),
##                                         ZmumuCandidatesCollection = cms.InputTag("ZmumuCandidatesMCJetHT")
## )

analyzeLumiMC = cms.EDAnalyzer("LumiAnalyzer",
                               genEventInfo = cms.InputTag("generator")
)

analyzeMuMuGen = cms.EDAnalyzer("MuMuGenAnalyzer",
                                GenParticleCollection = cms.InputTag("prunedGenParticles"),
                                GenJetCollection = cms.InputTag("slimmedGenJets"),
                                genEventInfo = cms.InputTag("generator"),
                                pileupSummaryInfo = cms.InputTag("slimmedAddPileupInfo"),
                                puDataFileName = cms.string("PileupHistogram-goldenJSON-13tev-2016-69200ub.root"),
                                puMCFileName = cms.string("PileupMC.root")
)

analyzeMuMuGenLooseIso = cms.EDAnalyzer("MuMuGenAnalyzer",
                                GenParticleCollection = cms.InputTag("prunedGenParticles"),
                                GenJetCollection = cms.InputTag("slimmedGenJets"),
                                genEventInfo = cms.InputTag("generator"),
                                pileupSummaryInfo = cms.InputTag("slimmedAddPileupInfo"),
                                puDataFileName = cms.string("PileupHistogram-goldenJSON-13tev-2016-69200ub.root"),
                                puMCFileName = cms.string("PileupMC.root")
)

analyzeMuMuGenLooseIsoDB = cms.EDAnalyzer("MuMuGenAnalyzer",
                                GenParticleCollection = cms.InputTag("prunedGenParticles"),
                                GenJetCollection = cms.InputTag("slimmedGenJets"),
                                genEventInfo = cms.InputTag("generator"),
                                pileupSummaryInfo = cms.InputTag("slimmedAddPileupInfo"),
                                puDataFileName = cms.string("PileupHistogram-goldenJSON-13tev-2016-69200ub.root"),
                                puMCFileName = cms.string("PileupMC.root")
)



## Sequences
makePatMuonsZmumuSelection = cms.Sequence(
    doubleMuonHLTTrigger
    + patMuonsAfterKinCutsDoubleMu
    + patMuonsAfterLooseIDDoubleMu
    + ZmumuCandidates
    + ZmumuCandidatesFilter
    + selectedMuonsForEmbedding
)

studyPatMuonsZmumuSelectionDoubleMu = cms.Sequence(
    doubleMuonHLTTrigger
    + patMuonsAfterKinCutsDoubleMu
    + patMuonsAfterLooseIDDoubleMu
    + ZmumuCandidates
    + ZmumuCandidatesFilter
    + analyzeMuonsForEmbedding
    + patMuonsAfterLooseIDLooseIsoDoubleMu
    + ZmumuCandidatesLooseIso
    + ZmumuCandidatesFilterLooseIso
    + analyzeMuonsForEmbeddingLooseIso
    + patMuonsAfterLooseIDLooseIsoDBDoubleMu
    + ZmumuCandidatesLooseIsoDB
    + ZmumuCandidatesFilterLooseIsoDB
    + analyzeMuonsForEmbeddingLooseIsoDB
)

studyPatMuonsZmumuSelectionSingleMu = cms.Sequence(
    singleMuonHLTTrigger
    + patMuonsAfterKinCutsSingleMu
    + patMuonsAfterLooseIDSingleMu
    + ZmumuCandidates
    + ZmumuCandidatesFilter
    + analyzeMuonsForEmbedding
)

studyPatMuonsZmumuSelectionJetHT = cms.Sequence(
    jetHTHLTTrigger
    + patMuonsAfterKinCutsJetHT
    + patMuonsAfterLooseIDJetHT
    + ZmumuCandidates
    + ZmumuCandidatesFilter
    + analyzeMuonsForEmbedding
)

studyPatMuonsZmumuSelectionMC = cms.Sequence(
    analyzeLumiMC 
    + doubleMuonHLTTrigger
    + patMuonsAfterKinCutsDoubleMu
    + patMuonsAfterLooseIDDoubleMu
    + ZmumuCandidates
    + ZmumuCandidatesFilter
    + analyzeMuonsForEmbedding
    + analyzeMuMuGen
    + patMuonsAfterLooseIDLooseIsoDoubleMu
    + ZmumuCandidatesLooseIso
    + ZmumuCandidatesFilterLooseIso
    + analyzeMuonsForEmbeddingLooseIso
    + analyzeMuMuGenLooseIso
    + patMuonsAfterLooseIDLooseIsoDBDoubleMu
    + ZmumuCandidatesLooseIsoDB
    + ZmumuCandidatesFilterLooseIsoDB
    + analyzeMuonsForEmbeddingLooseIsoDB
    + analyzeMuMuGenLooseIsoDB
)

## studyPatMuonsZmumuSelectionJetHT = cms.Sequence(
##     jetHTHLTTrigger
##     + patMuonsAfterKinCutsBoosted
##     + patMuonsAfterID
##     + patJetsAfterKinCuts
##     + patJetsAfterLooseID
##     + HighPtJetFilter
##     + ZmumuCandidates
##     + ZmumuCandidatesFilter
##     + analyzeMuonsForEmbedding
## )
