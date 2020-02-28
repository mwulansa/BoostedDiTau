import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')


options.maxEvents = -1
options.register('skipEvents', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Events to skip")
options.register('reportEvery', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Report every")
options.register('isMC', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Sample is MC")
options.register('doSlimming', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Output content is reduced")
options.register('isStandard', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Output content is reduced")
options.register('isReMINIAOD', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Output content is reduced")
options.register('numThreads', 8, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Set number of threads")

#inputFiles = 'root://cmseos.fnal.gov//eos/uscms/store/user/zhangj/events/ALP/RunIISummer17DR94Premix/ALP_m10_w1_htjmin400_RunIISummer17DR94Premix_AODSIM_100.root'

inputFiles = 'file:ALP_m10_w1_htjmin400_RunIISummer17DR94Premix_MINIAODSIM_Cleaned_100.root'

outputFile = 'ALP_m10_w1_htjmin400_RunIISummer17DR94Premix_reMINIAODSIM_Cleaned_100.root'

options.parseArguments()


#########################
### Main MINIAOD Path ###
#########################
from Configuration.StandardSequences.Eras import eras

process = cms.Process('reMINIAOD',eras.Run2_2017,eras.run2_miniAOD_94XFall17)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
if options.isMC:
    process.load('Configuration.StandardSequences.PATMC_cff')
else:
    process.load('Configuration.StandardSequences.PAT_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(inputFiles),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents = cms.untracked.uint32(options.skipEvents),
)

# Other statements
envvar = 'mcgt' if options.isMC else 'datagt'
from Configuration.AlCa.GlobalTag import GlobalTag
GT = {'mcgt': '94X_mc2017_realistic_v17', 'datagt': '94X_dataRun2_v11',}
process.GlobalTag = GlobalTag(process.GlobalTag, GT[envvar], '')

# New tauID
from BoostedDiTau.Skimmer.runTauIdMVA import *
na = TauIDEmbedder(process, cms, # pass tour process object
                   debug=True,
                   toKeep = ["2017v2"] # pick the one you need: ["2017v1", "2017v2", "newDM2017v2", "dR0p32017v2", "2016v1", "newDM2016v1"]
)
na.runTauID()

process.p = cms.Path(process.rerunMvaIsolationSequence
                     * process.NewTauIDsEmbedded ) # *getattr(process, "NewTauIDsEmbedded")


process.MINIAODoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('MINIAODSIM' if options.isMC else 'MINIAOD'),
        filterName = cms.untracked.string('')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    eventAutoFlushCompressedSize = cms.untracked.int32(-900),
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string(outputFile),
    outputCommands = process.MINIAODSIMEventContent.outputCommands,
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
process.MINIAODoutput_step = cms.EndPath(process.MINIAODoutput)

process.schedule = cms.Schedule(process.p, process.endjob_step, process.MINIAODoutput_step)


# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC, miniAOD_customizeAllData

#call to customisation function miniAOD_customizeAllMC imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
if options.isMC:
    process = miniAOD_customizeAllMC(process)
else:
    process = miniAOD_customizeAllData(process)


if not options.isStandard:
    from BoostedDiTau.Skimmer.customizeSkims4TCP import addCustomization
    addCustomization(process,options)
