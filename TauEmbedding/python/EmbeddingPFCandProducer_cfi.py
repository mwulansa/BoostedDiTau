import FWCore.ParameterSet.Config as cms


embeddedPFCandidates = cms.EDProducer("EmbeddingPFCandProducer",
                                      muonsSrc = cms.InputTag("selectedMuonsForEmbedding","","SELECT"),
                                      packedCandSrc = cms.InputTag("packedPFCandidates","","DQM"),
                                      packedCandSrcEmbedding = cms.InputTag("packedPFCandidates","","SIMembedding"),
                                      lostTrackSrc = cms.InputTag("lostTracks","","DQM"),
                                      lostTrackSrcEmbedding = cms.InputTag("lostTracks","","SIMembedding")
)

makeEmbeddedPFCandidates = cms.Sequence(embeddedPFCandidates)
