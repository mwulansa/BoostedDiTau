import FWCore.ParameterSet.Config as cms


externalLHEProducer = cms.EDProducer("EmbeddingLHEProducer",
    src = cms.InputTag("selectedMuonsForEmbedding","","SELECT"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices", "", ""),
    particleToEmbed = cms.int32(15),
    rotate180 = cms.bool(False),
    mirror = cms.bool(False),
    initialRecoCorrection = cms.bool(True),
    studyFSRmode = cms.untracked.string("reco")
)

makeexternalLHEProducer = cms.Sequence( externalLHEProducer)
