import FWCore.ParameterSet.Config as cms

flashggVertexMapUnique = cms.EDProducer('FlashggDzVertexMapProducer',
                                        PFCandidatesTag=cms.InputTag('packedPFCandidates'),
                                        VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                        MaxAllowedDz=cms.double(0.2),
                                        UseEachTrackOnce=cms.bool(True)
                                        )

flashggVertexMapUniqueZMuMu = cms.EDProducer('FlashggDzVertexMapProducer',
                                             PFCandidatesTag=cms.InputTag('packedPFCandidates'),
                                             VertexTag=cms.InputTag('offlinePrimaryVerticesWithMu'),
                                             MaxAllowedDz=cms.double(0.2),
                                             UseEachTrackOnce=cms.bool(True)
                                             )

flashggVertexMapUniqueZMuMuNoMu = cms.EDProducer('FlashggDzVertexMapProducer',
                                                 PFCandidatesTag=cms.InputTag('packedPFCandidates'),
                                                 VertexTag=cms.InputTag('offlinePrimaryVerticesNoMu'),
                                                 MaxAllowedDz=cms.double(0.2),
                                                 UseEachTrackOnce=cms.bool(True)
                                                 )

flashggVertexMapNonUnique = cms.EDProducer('FlashggDzVertexMapProducer',
                                           PFCandidatesTag=cms.InputTag('packedPFCandidates'),
                                           VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                           MaxAllowedDz=cms.double(0.2), 
                                           UseEachTrackOnce=cms.bool(False)
                                           )

flashggVertexMapNonUniqueZMuMu = cms.EDProducer('FlashggDzVertexMapProducer',
                                                PFCandidatesTag=cms.InputTag('packedPFCandidates'),
                                                VertexTag=cms.InputTag('offlinePrimaryVerticesWithMu'),
                                                MaxAllowedDz=cms.double(0.2), 
                                                UseEachTrackOnce=cms.bool(False)
                                                )

flashggVertexMapNonUniqueZMuMuNoMu = cms.EDProducer('FlashggDzVertexMapProducer',
                                                    PFCandidatesTag=cms.InputTag('packedPFCandidates'),
                                                    VertexTag=cms.InputTag('offlinePrimaryVerticesNoMu'),
                                                    MaxAllowedDz=cms.double(0.2), 
                                                    UseEachTrackOnce=cms.bool(False)
                                                    )
