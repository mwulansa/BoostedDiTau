import FWCore.ParameterSet.Config as cms
import sys
from Configuration.StandardSequences.Eras import eras

#inputFiles = "file:out_reco_pat.root"
#inputFiles = "root://xrootd.unl.edu//store/user/nbower/Events/TCP_m_30_w_1_htj_400toInf_slc6_amd64_gcc630_MINIAOD/TCP_m_30_w_1_htj_400toInf_slc6_amd64_gcc630_MINIAOD_97.root"
#outputFile = "embedded.root"
#maxEvents = -1


#prefix = "root://cmseos.fnal.gov//eos/uscms/store/group/lpcsusyhiggs/Jingyu/embedding/YMuMu/UL2017/"
prefix = "file:"

inputFile = prefix + sys.argv[2]
#inputFile = ''
outputFile = sys.argv[2].replace('simulated', 'embed')
#outputFile = ''
maxEvents = -1

process = cms.Process('Embed',eras.Run2_2017)

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

#process.load('BoostedDiTau.TauEmbedding.SelectingProcedure_cff')
#process.selection = cms.Path(process.makePatMuonsZmumuSelection)

process.MINIAODSIMoutput.outputCommands += [
    #'drop *_*_*_Embed',
    'keep *_*_*_SELECT',
    'keep *_*_*_LHEembedding',
    'keep *_*_*_SIMembedding',
    'keep *_selectedPatTausEmbed_*_*',
    'keep *_embeddedPFCandidates_*_*',
    'keep *_slimmedMuonsEmbedded_*_*',
    'keep *_selectedPatJetsAK4PFchs__*',
    'keep *_slimmedMETsTEST_*_*',
    'keep *_ak4PFJetsPATEmbed_*_*',
    'keep *_ak4PFJets_*_*',
    'keep *_patTausEmbed_*_*',
    'keep *_hpsPFTauProducerSansRefsEmbed_*_*',
    'keep *_combinatoricRecoTausEmbed_*_*',
    'keep *_ak4PFJetsRecoTauChargedHadronsEmbed_*_*',
    'keep *_ak4PFJetsLegacyHPSPiZerosEmbed_*_*',
    'keep *_patJetsPATEmbed_*_*',
    'keep *_ca8PFJetsCHSprunedForBoostedTausPATBoostedEmbed_*_*',
    'keep *_selectedPatTausNoNewIDsBoostedEmbed_*_*',
    'keep *_selectedPatTausBoostedEmbed_*_*'
    ]

process.endjob_step = cms.EndPath(process.endOfProcess)
process.MINIAODSIMoutput_step = cms.EndPath(process.MINIAODSIMoutput)

#process.load('BoostedDiTau.TauEmbedding.EmbeddingLHEProducer_cfi')
#process.lheproduction = cms.Path(process.makeexternalLHEProducer)

process.load('BoostedDiTau.TauEmbedding.EmbeddingPFCandProducer_cfi')
process.p= cms.Path(process.makeEmbeddedPFCandidates)

#process.p = cms.Path(process.makePatMuonsZmumuSelection + process.makeexternalLHEProducer)
#process.p = cms.Path(process.makePatMuonsZmumuSelection)

#print(process.embed)
#process.p = cms.Path(process.embed)

#process.schedule = cms.Schedule(process.p, process.endjob_step, process.MINIAODoutput_step)

#process.schedule.associate(process.patTask)
#from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
#associatePatAlgosToolsTask(process)


from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
jetToolbox( process, 'ak4', 'dummy', 'out', PUMethod='chs', dataTier="miniAOD", bTagDiscriminators=["pfDeepCSVJetTags:probb","pfDeepCSVJetTags:probbb"], newPFCollection=True, nameNewPFCollection="embeddedPFCandidates")

process.ak4PFJets.src = cms.InputTag("embeddedPFCandidates","packedPFCandidatesEmbedded")
#process.ak4PFJets.srcPVs = cms.InputTag("")

process.pfInclusiveSecondaryVertexFinderTagInfosAK4PFchs.extSVCollection = cms.InputTag("slimmedSecondaryVertices","", "PAT")

#process.candidateVertexArbitratorCvsL.tracks = cms.InputTag("embeddedPFCandidates","lostTracksEmbedded")
#process.inclusiveCandidateVertexFinderCvsL.tracks = cms.InputTag("embeddedPFCandidates","lostTracksEmbedded")
process.pfImpactParameterTagInfosAK4PFchs.candidates = cms.InputTag("embeddedPFCandidates","packedPFCandidatesEmbedded")

#####
from RecoTauTag.Configuration.tools.adaptToRunAtMiniAOD import adaptToRunAtMiniAODEmbedding

postfix = 'Embed'
runBoosted = True

#####
tauAtMiniTools = adaptToRunAtMiniAODEmbedding(process,runBoosted,postfix=postfix)
tauAtMiniTools.addTauReReco()
tauAtMiniTools.adaptTauToMiniAODReReco()
process.p1 = cms.Path(
        getattr(process,("miniAODTausSequence"+postfix if not runBoosted else "miniAODTausSequenceBoosted"+postfix))
)


#######
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

runMetCorAndUncFromMiniAOD(process,
                           isData=False,
                           pfCandColl=cms.InputTag("embeddedPFCandidates","packedPFCandidatesEmbedded"),
                           recoMetFromPFCs=True,
                           CHS = True, #This is an important step and determines what type of jets to be reclustered
                           reclusterJets = True,
                           postfix="TEST"
                           )

process.packedPrimaryVertexAssociationJME.jets = cms.InputTag("selectedPatJetsAK4PFchs","","Embed")
process.packedPrimaryVertexAssociationJME.particles = cms.InputTag("embeddedPFCandidates","packedPFCandidatesEmbedded")
process.packedPrimaryVertexAssociationJME.vertices = cms.InputTag("offlineSlimmedPrimaryVertices", "", "PAT")

process.pfCHS.candidates = cms.InputTag("embeddedPFCandidates","packedPFCandidatesEmbedded")
#process.pfCHS.vertexAssociation = cms.InputTag("packedPrimaryVertexAssociationJME","","Embed")

process.basicJetsForMetTEST.jetCorrLabelRes = cms.InputTag("L3Absolute")
#process.basicJetsForMetNoHFTEST.jetCorrLabelRes = cms.InputTag("L3Absolute")
#process.basicJetsForMetPuppiTEST.jetCorrLabelRes = cms.InputTag("L3Absolute")

process.patPFMetT1T2CorrTEST.jetCorrLabelRes = cms.InputTag("L3Absolute")
#process.patPFMetT1T2CorrNoHFTEST.jetCorrLabelRes = cms.InputTag("L3Absolute")
#process.patPFMetT1T2CorrPuppiTEST.jetCorrLabelRes = cms.InputTag("L3Absolute")

process.patPFMetT1T2SmearCorrTEST.jetCorrLabelRes = cms.InputTag("L3Absolute")
#process.patPFMetT1T2SmearCorrNoHFTEST.jetCorrLabelRes = cms.InputTag("L3Absolute")
#process.patPFMetT1T2SmearCorrPuppiTEST.jetCorrLabelRes = cms.InputTag("L3Absolute")


#process.ak4PFJetsLegacyHPSPiZeros.builders.verbosity = cms.int32(1)

#process.pfPileUpIsoPFBRECO.PFCandidates = cms.InputTag("embeddedPFCandidates")
#process.pfPileUpIsoPFBRECO.Enable = cms.bool(False)
#process.pfNoPileUpIsoPFBRECO.enable = cms.bool(False)

#process.patJetPartonsLegacy.src = cms.InputTag("prunedGenParticles")
#process.genParticlesForJetsNoNuTMP.src = cms.InputTag("prunedGenParticles")

# Customisation from command line
#process.RandomNumberGeneratorService.externalLHEProducer.initialSeed=int(1552064665%100)


process.deepTau2017v2p1MiniAODTausBoostedEmbed.electrons = cms.InputTag("slimmedElectrons","","PAT")
process.deepTau2017v2p1MiniAODTausBoostedEmbed.muons = cms.InputTag('embeddedPFCandidates','slimmedMuonsEmbedded')
process.deepTau2017v2p1MiniAODTausBoostedEmbed.pfcands = cms.InputTag('embeddedPFCandidates','packedPFCandidatesEmbedded')
process.deepTau2017v2p1MiniAODTausBoostedEmbed.vertices = cms.InputTag('embeddedPFCandidates', 'offlineSlimmedPrimaryVerticesEmbedded')

process.hpsPFTauPrimaryVertexProducerBoostedEmbed.lostCandidatesTag = cms.InputTag('embeddedPFCandidates','lostTracksEmbedded')

process.embeddedPFCandidates.muonsSrc = cms.InputTag("selectedBoostedMuonsForEmbedding","","SELECT")
process.embeddedPFCandidates.muonsEmbedSrc = cms.InputTag("slimmedMuons","","PAT")
process.embeddedPFCandidates.electronsSrc = cms.InputTag("slimmedElectrons","","PAT")
process.embeddedPFCandidates.packedCandSrc = cms.InputTag("packedPFCandidates","","PAT")
process.embeddedPFCandidates.lostTrackSrc = cms.InputTag("lostTracks","","PAT")
process.embeddedPFCandidates.offlineSlimmedPrimaryVerticesSrc = cms.InputTag("offlineSlimmedPrimaryVertices","", "PAT")
