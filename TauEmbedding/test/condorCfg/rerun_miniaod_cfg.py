import FWCore.ParameterSet.Config as cms
import sys
from Configuration.StandardSequences.Eras import eras

#inputFiles = "file:out_reco_pat.root"
#inputFiles = "root://xrootd.unl.edu//store/user/nbower/Events/TCP_m_30_w_1_htj_400toInf_slc6_amd64_gcc630_MINIAOD/TCP_m_30_w_1_htj_400toInf_slc6_amd64_gcc630_MINIAOD_97.root"
#outputFile = "embedded.root"
#maxEvents = -1


prefix = "root://cmseos.fnal.gov//store/user/zhangj/events/UpsilonTauTau/UL2017/"
#prefix = "file:"

#inputFile = prefix + sys.argv[2]
inputFile = ''
#outputFile = 'reminiaod_'+sys.argv[2]
outputFile = ''
maxEvents = -1

process = cms.Process('ReMINIAOD',eras.Run2_2017)

# import of standard configurations
#process.load('Configuration.StandardSequences.Services_cff')
#process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.load('Configuration.StandardSequences.PATMC_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(maxEvents)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(inputFile),
    secondaryFileNames = cms.untracked.vstring(),
)

#from Configuration.AlCa.GlobalTag import GlobalTag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_realistic', '')


process.MINIAODSIMoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('MINIAODSIM'),
        filterName = cms.untracked.string('')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    eventAutoFlushCompressedSize = cms.untracked.int32(-900),
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string(outputFile),
    outputCommands = process.MINIAODSIMEventContent.outputCommands,
    SelectEvents=cms.untracked.PSet(
        SelectEvents=cms.vstring('p')
    ),
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

process.endjob_step = cms.EndPath(process.endOfProcess)
process.MINIAODoutput_step = cms.EndPath(process.MINIAODSIMoutput)

#process.schedule = cms.Schedule(process.endjob_step, process.MINIAODoutput_step)

process.MINIAODSIMoutput.outputCommands += [
    'keep *_ca8PFJetsCHSprunedForBoostedTausPATBoostedReMINIAOD_*_*',
    'keep *_selectedPatTausNoNewIDsBoostedReMINIAOD_*_*',
    'keep *_selectedPatTausBoostedReMINIAOD_*_*', 
    'keep *_BoostedReMINIAOD_*_*'
    ]

#process.schedule.associate(process.patTask)
#from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
#associatePatAlgosToolsTask(process)


#####
from RecoTauTag.Configuration.tools.adaptToRunAtMiniAOD import adaptToRunAtMiniAOD

postfix = 'ReMINIAOD'
runBoosted = True

#####
tauAtMiniTools = adaptToRunAtMiniAOD(process,runBoosted,postfix=postfix)
tauAtMiniTools.addTauReReco()
tauAtMiniTools.adaptTauToMiniAODReReco()
process.p = cms.Path(
        getattr(process,("miniAODTausSequence"+postfix if not runBoosted else "miniAODTausSequenceBoosted"+postfix))
)


