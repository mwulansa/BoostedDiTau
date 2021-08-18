import os, sys, random
import numpy as np

configDir="./configs/"

#os.mkdir(configDir)

storagePrefix1='root://cmseos.fnal.gov//eos/uscms/store/user/nbower/Events/ALP/'
storagePrefix2='root://cmseos.fnal.gov//eos/uscms/store/user/zhangj/events/ALP/RunIISummer19UL17RECO/'
#storagePrefix2 = 'root://cmseos.fnal.gov//eos/uscms/store/user/zhangj/events/UpsilonTauTau/RunIISummer19UL17RECO/'

#masses=[10, 30, 50]
#masses=[10]
masses = [30, 50]
jobs=np.linspace(100,1,100)

for mass in masses:
    for job in jobs:
        configFile = "ALP_m"+str(mass)+"_w1_htjmin400_RunIISummer19UL17RECO_MINIAODSIM_"+str(int(job))+".py"
        inputFile = "TCP_m"+str(mass)+"_w1_htjmin400_RunIISummer19UL17RECO_AODSIM_"+str(int(job))+".root"
        #configFile = "UpsilonToTauTau_pthatmin400_RunIISummer19UL17RECO_MINIAODSIM_"+str(int(job))+".py"
        #inputFile = "UpsilonToTauTau_pthatmin400_RunIISummer19UL17RECO_AODSIM_"+str(int(job))+".root"
        print configDir+configFile
        cfg=open(configDir+configFile, "w")
        cfg.writelines("""
import FWCore.ParameterSet.Config as cms


from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')

inputFiles = '"""+storagePrefix2+inputFile+"""'

options.maxEvents = -1
options.register('skipEvents', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Events to skip")
options.register('reportEvery', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Report every")
options.register('isMC', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Sample is MC")
options.register('doSlimming', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Output content is reduced")
options.register('isStandard', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Output content is reduced")
options.register('isReMINIAOD', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Output content is reduced")
options.register('numThreads', 8, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Set number of threads")

if options.isStandard:
    options.doSlimming = 0
    outputFile = '"""+storagePrefix2+configFile.replace("MINIAODSIM", "MINIAODSIM_Standard").replace(".py",".root")+"""'
elif options.doSlimming:
    outputFile = '"""+storagePrefix2+configFile.replace("MINIAODSIM", "MINIAODSIM_Slimmed").replace(".py",".root")+"""'
else:
    outputFile = '"""+storagePrefix2+configFile.replace("MINIAODSIM", "MINIAODSIM_Cleaned").replace(".py",".root")+"""'

options.parseArguments()

#########################
### Main MINIAOD Path ###
#########################
from Configuration.Eras.Era_Run2_2017_cff import Run2_2017

process = cms.Process('PAT',Run2_2017)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.PATMC_cff')
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

process.options = cms.untracked.PSet(
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step2 nevts:4800'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition
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
    overrideBranchesSplitLevel = cms.untracked.VPSet(
        cms.untracked.PSet(
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
        )
    ),
    overrideInputFileSplitLevels = cms.untracked.bool(True),
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mc2017_realistic_v6', '')

#patAlgosToolsTask.add(process.selectedPATTausElectronCleaned)

# Path and EndPath definitions
#process.tauId = cms.Path(
#    process.rerunMvaIsolationSequence *
#    getattr(process,updatedTauName)
#)

# Path and EndPath definitions
process.Flag_trackingFailureFilter = cms.Path(process.goodVertices+process.trackingFailureFilter)
process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)
process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)
process.Flag_trkPOGFilters = cms.Path(process.trkPOGFilters)
process.Flag_HcalStripHaloFilter = cms.Path(process.HcalStripHaloFilter)
process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(~process.logErrorTooManyClusters)
process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
process.Flag_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)
process.Flag_globalSuperTightHalo2016Filter = cms.Path(process.globalSuperTightHalo2016Filter)
process.Flag_eeBadScFilter = cms.Path(process.eeBadScFilter)
process.Flag_METFilters = cms.Path(process.metFilters)
process.Flag_chargedHadronTrackResolutionFilter = cms.Path(process.chargedHadronTrackResolutionFilter)
process.Flag_globalTightHalo2016Filter = cms.Path(process.globalTightHalo2016Filter)
process.Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(process.CSCTightHaloTrkMuUnvetoFilter)
process.Flag_HBHENoiseIsoFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseIsoFilter)
process.Flag_BadChargedCandidateSummer16Filter = cms.Path(process.BadChargedCandidateSummer16Filter)
process.Flag_hcalLaserEventFilter = cms.Path(process.hcalLaserEventFilter)
process.Flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter)
process.Flag_ecalBadCalibFilter = cms.Path(process.ecalBadCalibFilter)
process.Flag_HBHENoiseFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseFilter)
process.Flag_trkPOG_toomanystripclus53X = cms.Path(~process.toomanystripclus53X)
process.Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path(process.EcalDeadCellBoundaryEnergyFilter)
process.Flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter)
process.Flag_trkPOG_manystripclus53X = cms.Path(~process.manystripclus53X)
process.Flag_BadPFMuonSummer16Filter = cms.Path(process.BadPFMuonSummer16Filter)
process.Flag_muonBadTrackFilter = cms.Path(process.muonBadTrackFilter)
process.Flag_CSCTightHalo2015Filter = cms.Path(process.CSCTightHalo2015Filter)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.MINIAODSIMoutput_step = cms.EndPath(process.MINIAODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.Flag_HBHENoiseFilter,process.Flag_HBHENoiseIsoFilter,process.Flag_CSCTightHaloFilter,process.Flag_CSCTightHaloTrkMuUnvetoFilter,process.Flag_CSCTightHalo2015Filter,process.Flag_globalTightHalo2016Filter,process.Flag_globalSuperTightHalo2016Filter,process.Flag_HcalStripHaloFilter,process.Flag_hcalLaserEventFilter,process.Flag_EcalDeadCellTriggerPrimitiveFilter,process.Flag_EcalDeadCellBoundaryEnergyFilter,process.Flag_ecalBadCalibFilter,process.Flag_goodVertices,process.Flag_eeBadScFilter,process.Flag_ecalLaserCorrFilter,process.Flag_trkPOGFilters,process.Flag_chargedHadronTrackResolutionFilter,process.Flag_muonBadTrackFilter,process.Flag_BadChargedCandidateFilter,process.Flag_BadPFMuonFilter,process.Flag_BadChargedCandidateSummer16Filter,process.Flag_BadPFMuonSummer16Filter,process.Flag_trkPOG_manystripclus53X,process.Flag_trkPOG_toomanystripclus53X,process.Flag_trkPOG_logErrorTooManyClusters,process.Flag_METFilters,process.endjob_step,process.MINIAODSIMoutput_step)
process.schedule.associate(process.patTask)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(options.numThreads)
process.options.numberOfStreams=cms.untracked.uint32(0)
process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)

#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC 

#call to customisation function miniAOD_customizeAllMC imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
process = miniAOD_customizeAllMC(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

# Add the customization
if not options.isStandard:
    from BoostedDiTau.Skimmer.customizeSkims4TCP import addCustomization
    addCustomization(process,options)
## fix things
process.patTrigger.processName = cms.string('RECO')
process.slimmedPatTrigger.triggerResults = cms.InputTag("TriggerResults","","RECO")

# add tauID
updatedTauName = "slimmedTausNewID"
#import RecoTauTag.RecoTau.tools.runTauIdMVA as tauIdConfig
import BoostedDiTau.Skimmer.runTauIdMVA as tauIdConfig
tauIdEmbedder = tauIdConfig.TauIDEmbedder(process, cms, debug = True,
                                          updatedTauName = updatedTauName,
                                          tauSrc = 'slimmedTaus',
                                          toKeep = ["2017v2",
                                                    "deepTau2017v2p1", #deepTau TauIDs
                                          ])
tauIdEmbedder.runTauID()

import PhysicsTools.PatAlgos.tools.helpers as configtools
patAlgosToolsTask = configtools.getPatAlgosToolsTask(process)

process.tauId = cms.Task(
    process.rerunMvaIsolationTask,
    process.slimmedTausNewID,
)

patAlgosToolsTask.add(process.tauId)

process.MINIAODSIMEventContent.outputCommands += [
    'keep *_slimmedTausNewID_*_*',
]
        """)
