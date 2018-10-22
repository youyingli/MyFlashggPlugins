import FWCore.ParameterSet.Config as cms
from MyFlashggPlugins.DiphotonVertex.flashggTkVtxMapValidation_cfi import flashggVertexMapUniqueZMuMu,flashggVertexMapUniqueZMuMuNoMu,flashggVertexMapNonUniqueZMuMu,flashggVertexMapNonUniqueZMuMuNoMu

flashggZMuMuValidationSequence = cms.Sequence(flashggVertexMapUniqueZMuMu + flashggVertexMapUniqueZMuMuNoMu
                                            + flashggVertexMapNonUniqueZMuMu + flashggVertexMapNonUniqueZMuMuNoMu
                                             )
