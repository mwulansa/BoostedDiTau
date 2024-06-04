import FWCore.ParameterSet.Config as cms


embeddedPFCandidates = cms.EDProducer("EmbeddingPFCandProducer",
                                      muonsSrc = cms.InputTag("selectedMuonsForEmbedding","","SELECT"),
                                      muonsEmbedSrc = cms.InputTag("slimmedMuons","","DQM"),
                                      muonsSrcEmbedding = cms.InputTag("slimmedMuons","","SIMembedding"),
                                      electronsSrc = cms.InputTag("slimmedElectrons","","DQM"),
                                      electronsSrcEmbedding = cms.InputTag("slimmedElectrons","","SIMembedding"),
                                      packedCandSrc = cms.InputTag("packedPFCandidates","","DQM"),
                                      packedCandSrcEmbedding = cms.InputTag("packedPFCandidates","","SIMembedding"),
                                      lostTrackSrc = cms.InputTag("lostTracks","","DQM"),
                                      lostTrackSrcEmbedding = cms.InputTag("lostTracks","","SIMembedding"),
                                      offlineSlimmedPrimaryVerticesSrc = cms.InputTag("offlineSlimmedPrimaryVertices","", "DQM"),
                                      offlineSlimmedPrimaryVerticesSrcEmbedding = cms.InputTag("offlineSlimmedPrimaryVertices","","SIMembedding")
)

makeEmbeddedPFCandidates = cms.Sequence(embeddedPFCandidates)
