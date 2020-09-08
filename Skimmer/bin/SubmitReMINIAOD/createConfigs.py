import os, sys, random
import numpy as np

configDir="./configs/"

#os.mkdir(configDir)

inputFileListName = sys.argv[1]
sample = sys.argv[2]
mass = sys.argv[3]

inputFileList = inputFileListName

OutputDir = 'root://cmseos.fnal.gov//store/user/mwulansa/DIS/TCPAnalysis/Backgrounds/RunIIFall17DR94Premix/'

inputSampleName = inputFileListName.replace("filelists/"+sample+"/"+mass+"/", "")

print configDir+inputSampleName

cfg=open (configDir+inputSampleName.replace(".txt", ".py"), "w")
cfg.writelines("""


from PhysicsTools.PatAlgos.tools.helpers import cloneProcessingSnippet
from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')

options.maxEvents = 10000
options.register('skipEvents', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Events to skip")
options.register('reportEvery', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Report every")
options.register('isMC', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Sample is MC")
options.register('doSlimming', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Output content is reduced")
options.register('isStandard', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Output content is reduced")
options.register('isReMINIAOD', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Output content is reduced")
options.register('numThreads', 8, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Set number of threads")

if options.isStandard:
    options.doSlimming = 0
    outputFile = '"""+OutputDir+"""standard_"""+inputSampleName.replace(".txt",".root"\
)+"""'
elif options.doSlimming:
    outputFile = '"""+OutputDir+"""slimmed_"""+inputSampleName.replace(".txt",".root")\
+"""'
else:
    outputFile = '"""+OutputDir+inputSampleName.replace(".txt",".root")+"""'

options.parseArguments()


#########################
### Main MINIAOD Path ###
#########################
from Configuration.StandardSequences.Eras import eras

process = cms.Process('PAT',eras.Run2_2017,eras.run2_miniAOD_94XFall17)

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
    fileNames = cms.untracked.vstring(

""")

cfg.close()

inputFileNames = open(inputFileList, 'r')
for inputFileName in inputFileNames:
#    print inputFileName
    cfg = open (configDir+inputSampleName.replace(".txt", ".py"), "a")
#    cfg = open (configDir+inputFileListName.replace(".txt", ".py"), "a")
#    cfg.writelines("""'file:"""+inputFileName.replace("\n","")+"""',\n""")
    cfg.writelines("""'"""+inputFileName.replace("\n","")+"""',\n""")
#    cfg.writelines("""'"""+inputFileName.replace("\n","")+"""',\n""")

cfg.close()
cfg = open (configDir+inputSampleName.replace(".txt", ".py"), "a")
#cfg = open (configDir+inputFileListName.replace(".txt", ".py"), "a")
cfg.writelines("""
),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents = cms.untracked.uint32(options.skipEvents),
)

# Other statements
envvar = 'mcgt' if options.isMC else 'datagt'
from Configuration.AlCa.GlobalTag import GlobalTag
GT = {'mcgt': '94X_mc2017_realistic_v17', 'datagt': '94X_dataRun2_v11',}
process.GlobalTag = GlobalTag(process.GlobalTag, GT[envvar], '')

# New TauID
from BoostedDiTau.Skimmer.runTauIdMVA import *
na = TauIDEmbedder(process, cms, # pass tour process object
                   debug=True,
                   toKeep = ["2017v2"] # pick the one you need: ["2017v1", "2017v2", "newDM2017v2", "dR0p32017v2", "2016v1", "newDM2016v1"]
)

#na.runTauID('NewTauIDsEmbedded')
na.runTauID()

process.rerunMvaIsolationSequenceMuonCleaned=cloneProcessingSnippet(process,process.rerunMvaIsolationSequence,'MuonCleaned', addToTask =False)
massSearchReplaceAnyInputTag(process.rerunMvaIsolationSequenceMuonCleaned,cms.InputTag('slimmedTaus'),cms.InputTag('slimmedTausMuonCleaned'))

process.rerunMvaIsolationSequenceElectronCleaned=cloneProcessingSnippet(process,process.rerunMvaIsolationSequence,'ElectronCleaned', addToTask =False)
massSearchReplaceAnyInputTag(process.rerunMvaIsolationSequenceElectronCleaned,cms.InputTag('slimmedTaus'),cms.InputTag('slimmedTausElectronCleaned'))

process.rerunMvaIsolationSequenceElectronCleanedFlat=cloneProcessingSnippet(process,process.rerunMvaIsolationSequence,'ElectronCleanedFlat', addToTask =False)
massSearchReplaceAnyInputTag(process.rerunMvaIsolationSequenceElectronCleanedFlat,cms.InputTag('slimmedTaus'),cms.InputTag('slimmedTausElectronCleanedFlat'))


#process.NewTauIDsEmbeddedMuonCleaned=cloneProcessingSnippet(process,process.NewTauIDsEmbedded,'MuonCleaned', addToTask =False)
#massSearchReplaceAnyInputTag(process.NewTauIDsEmbeddedMuonCleaned,cms.InputTag('slimmedTaus'),cms.InputTag('slimmedTausMuonCleaned'))

process.p = cms.Path(process.rerunMvaIsolationSequence
                     * process.NewTauIDsEmbedded
                     )

process.NewTauIDsEmbeddedMuonCleaned = process.NewTauIDsEmbedded.clone()
process.NewTauIDsEmbeddedMuonCleaned.src = cms.InputTag('slimmedTausMuonCleaned')

process.NewTauIDsEmbeddedElectronCleaned = process.NewTauIDsEmbedded.clone()
process.NewTauIDsEmbeddedElectronCleaned.src = cms.InputTag('slimmedTausElectronCleaned')

process.NewTauIDsEmbeddedElectronCleanedFlat = process.NewTauIDsEmbedded.clone()
process.NewTauIDsEmbeddedElectronCleanedFlat.src = cms.InputTag('slimmedTausElectronCleanedFlat')

process.NewTauIDsEmbeddedMuonCleaned.tauIDSources = cms.PSet(
        byIsolationMVArun2017v2DBoldDMwLTraw2017 = cms.InputTag("rerunDiscriminationByIsolationOldDMMVArun2017v2rawMuonCleaned"),
        byLooseIsolationMVArun2017v2DBoldDMwLT2017 = cms.InputTag("rerunDiscriminationByIsolationOldDMMVArun2017v2LooseMuonCleaned"),
        byMediumIsolationMVArun2017v2DBoldDMwLT2017 = cms.InputTag("rerunDiscriminationByIsolationOldDMMVArun2017v2MediumMuonCleaned"),
        byTightIsolationMVArun2017v2DBoldDMwLT2017 = cms.InputTag("rerunDiscriminationByIsolationOldDMMVArun2017v2TightMuonCleaned"),
        byVLooseIsolationMVArun2017v2DBoldDMwLT2017 = cms.InputTag("rerunDiscriminationByIsolationOldDMMVArun2017v2VLooseMuonCleaned"),
        byVTightIsolationMVArun2017v2DBoldDMwLT2017 = cms.InputTag("rerunDiscriminationByIsolationOldDMMVArun2017v2VTightMuonCleaned"),
        byVVLooseIsolationMVArun2017v2DBoldDMwLT2017 = cms.InputTag("rerunDiscriminationByIsolationOldDMMVArun2017v2VVLooseMuonCleaned"),
        byVVTightIsolationMVArun2017v2DBoldDMwLT2017 = cms.InputTag("rerunDiscriminationByIsolationOldDMMVArun2017v2VVTightMuonCleaned")
    )

process.NewTauIDsEmbeddedElectronCleaned.tauIDSources = cms.PSet(
        byIsolationMVArun2017v2DBoldDMwLTraw2017 = cms.InputTag("rerunDiscriminationByIsolationOldDMMVArun2017v2rawElectronCleaned"),
        byLooseIsolationMVArun2017v2DBoldDMwLT2017 = cms.InputTag("rerunDiscriminationByIsolationOldDMMVArun2017v2LooseElectronCleaned"),
        byMediumIsolationMVArun2017v2DBoldDMwLT2017 = cms.InputTag("rerunDiscriminationByIsolationOldDMMVArun2017v2MediumElectronCleaned"),
        byTightIsolationMVArun2017v2DBoldDMwLT2017 = cms.InputTag("rerunDiscriminationByIsolationOldDMMVArun2017v2TightElectronCleaned"),
        byVLooseIsolationMVArun2017v2DBoldDMwLT2017 = cms.InputTag("rerunDiscriminationByIsolationOldDMMVArun2017v2VLooseElectronCleaned"),
        byVTightIsolationMVArun2017v2DBoldDMwLT2017 = cms.InputTag("rerunDiscriminationByIsolationOldDMMVArun2017v2VTightElectronCleaned"),
        byVVLooseIsolationMVArun2017v2DBoldDMwLT2017 = cms.InputTag("rerunDiscriminationByIsolationOldDMMVArun2017v2VVLooseElectronCleaned"),
        byVVTightIsolationMVArun2017v2DBoldDMwLT2017 = cms.InputTag("rerunDiscriminationByIsolationOldDMMVArun2017v2VVTightElectronCleaned")
    )

process.NewTauIDsEmbeddedElectronCleanedFlat.tauIDSources = cms.PSet(
        byIsolationMVArun2017v2DBoldDMwLTraw2017 = cms.InputTag("rerunDiscriminationByIsolationOldDMMVArun2017v2rawElectronCleanedFlat"),
        byLooseIsolationMVArun2017v2DBoldDMwLT2017 = cms.InputTag("rerunDiscriminationByIsolationOldDMMVArun2017v2LooseElectronCleanedFlat"),
        byMediumIsolationMVArun2017v2DBoldDMwLT2017 = cms.InputTag("rerunDiscriminationByIsolationOldDMMVArun2017v2MediumElectronCleanedFlat"),
        byTightIsolationMVArun2017v2DBoldDMwLT2017 = cms.InputTag("rerunDiscriminationByIsolationOldDMMVArun2017v2TightElectronCleanedFlat"),
        byVLooseIsolationMVArun2017v2DBoldDMwLT2017 = cms.InputTag("rerunDiscriminationByIsolationOldDMMVArun2017v2VLooseElectronCleanedFlat"),
        byVTightIsolationMVArun2017v2DBoldDMwLT2017 = cms.InputTag("rerunDiscriminationByIsolationOldDMMVArun2017v2VTightElectronCleanedFlat"),
        byVVLooseIsolationMVArun2017v2DBoldDMwLT2017 = cms.InputTag("rerunDiscriminationByIsolationOldDMMVArun2017v2VVLooseElectronCleanedFlat"),
        byVVTightIsolationMVArun2017v2DBoldDMwLT2017 = cms.InputTag("rerunDiscriminationByIsolationOldDMMVArun2017v2VVTightElectronCleanedFlat")
    )


process.pm = cms.Path(process.rerunMvaIsolationSequenceMuonCleaned
                     * process.NewTauIDsEmbeddedMuonCleaned
                     )

process.pe = cms.Path(process.rerunMvaIsolationSequenceElectronCleaned
                     * process.NewTauIDsEmbeddedElectronCleaned
                     )

process.pef = cms.Path(process.rerunMvaIsolationSequenceElectronCleanedFlat
                     * process.NewTauIDsEmbeddedElectronCleanedFlat
                     )


#------------------------------------------------------------

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    #wantSummary = cms.untracked.bool(True),
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step1 nevts:4800'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)
#-------------------------------------------------------------

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

#--------------------------------------------------------

process.Flag_trackingFailureFilter = cms.Path(process.goodVertices+process.trackingFailureFilter)
process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)
process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)
process.Flag_trkPOGFilters = cms.Path(process.trkPOGFilters)
process.Flag_HcalStripHaloFilter = cms.Path(process.HcalStripHaloFilter)
process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(~process.logErrorTooManyClusters)
process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
process.Flag_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)
process.Flag_globalSuperTightHalo2016Filter = cms.Path(process.globalSuperTightHalo2016Filter)
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

#-------------------------------------------------------

process.endjob_step = cms.EndPath(process.endOfProcess)
process.MINIAODoutput_step = cms.EndPath(process.MINIAODoutput)

#---------------------------------------------------------------

process.schedule = cms.Schedule(process.p, process.pm, process.pe, process.pef, process.Flag_HBHENoiseFilter,process.Flag_HBHENoiseIsoFilter,process.Flag_CSCTightHaloFilter,process.Flag_CSCTightHaloTrkMuUnvetoFilter,process.Flag_CSCTightHalo2015Filter,process.Flag_globalTightHalo2016Filter,process.Flag_globalSuperTightHalo2016Filter,process.Flag_HcalStripHaloFilter,process.Flag_hcalLaserEventFilter,process.Flag_EcalDeadCellTriggerPrimitiveFilter,process.Flag_EcalDeadCellBoundaryEnergyFilter,process.Flag_ecalBadCalibFilter,process.Flag_goodVertices,process.Flag_eeBadScFilter,process.Flag_ecalLaserCorrFilter,process.Flag_trkPOGFilters,process.Flag_chargedHadronTrackResolutionFilter,process.Flag_muonBadTrackFilter,process.Flag_BadChargedCandidateFilter,process.Flag_BadPFMuonFilter,process.Flag_BadChargedCandidateSummer16Filter,process.Flag_BadPFMuonSummer16Filter,process.Flag_trkPOG_manystripclus53X,process.Flag_trkPOG_toomanystripclus53X,process.Flag_trkPOG_logErrorTooManyClusters,process.Flag_METFilters,process.endjob_step,process.MINIAODoutput_step)

process.schedule.associate(process.patTask)

#---------------------------------------------------------------

#process.schedule = cms.Schedule(process.p, process.pm, process.pe, process.endjob_step, process.MINIAODoutput_step)


# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC, miniAOD_customizeAllData

#--------------------------------
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

process.options.numberOfThreads=cms.untracked.uint32(options.numThreads)
process.options.numberOfStreams=cms.untracked.uint32(0)

from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)


#--------------------------------

#call to customisation function miniAOD_customizeAllMC imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
if options.isMC:
    process = miniAOD_customizeAllMC(process)
else:
    process = miniAOD_customizeAllData(process)


if not options.isStandard:
    from BoostedDiTau.Skimmer.customizeSkims4TCP import addCustomization
    addCustomization(process,options)

process.patTrigger.processName = cms.string('RECO')
process.slimmedPatTrigger.triggerResults = cms.InputTag("TriggerResults","","RECO")

dump_file = open('dump_config.py','w')
dump_file.write(process.dumpPython())

        """)
