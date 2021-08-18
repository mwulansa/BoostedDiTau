import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

inputFiles = "root://xrootd.unl.edu//store/data/Run2016B/DoubleMuon/MINIAOD/17Jul2018_ver2-v1/00000/0AB088EE-EA8A-E811-8636-0CC47A4C8E96.root"
outputFile = "outMuonSelection.root"
maxEvents = -1

process = cms.Process('reMINIAOD',eras.Run2_2017,eras.run2_miniAOD_94XFall17)

# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(maxEvents)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(inputFiles),
    secondaryFileNames = cms.untracked.vstring(),
)

process.MINIAODoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('MINIAOD'),
        filterName = cms.untracked.string('')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    eventAutoFlushCompressedSize = cms.untracked.int32(-900),
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string(outputFile),
    outputCommands = process.MINIAODEventContent.outputCommands,
    overrideBranchesSplitLevel = cms.untracked.VPSet(cms.untracked.PSet(
        branch = cms.untracked.string('patPackedCandidates_packedPFCandidates__*'),
        splitLevel = cms.untracked.int32(99)
    ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenParticles_prunedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patTriggerObjectStandAlones_slimmedPatTrigger__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patPackedGenParticles_packedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJets__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoVertexs_offlineSlimmedPrimaryVertices__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoCaloClusters_reducedEgamma_reducedESClusters_*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEBRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEERecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenJets_slimmedGenJets__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJetsPuppi__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedESRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        )),
    overrideInputFileSplitLevels = cms.untracked.bool(True),
    splitLevel = cms.untracked.int32(0)
)

process.load('BoostedDiTau.TauEmbedding.SelectingProcedure_cff')
process.p = cms.Path(process.makePatMuonsZmumuSelection)

process.MINIAODEventContent.outputCommands += [
    'keep *_patMuonsAfterID_*_*',
    #'keep *_ZmumuCandidates_*_*',
    'keep *_selectedMuonsForEmbedding_*_*',
    'keep *_selectedBoostedMuonsForEmbedding_*_*',
    ]

process.endjob_step = cms.EndPath(process.endOfProcess)
process.MINIAODoutput_step = cms.EndPath(process.MINIAODoutput)

process.schedule = cms.Schedule(process.p, process.endjob_step, process.MINIAODoutput_step)
