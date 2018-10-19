import FWCore.ParameterSet.Config as cms
from flashggDiphotonVtxPlugins.Validation.flashggTkVtxMapValidation_cfi import flashggVertexMapUnique,flashggVertexMapNonUnique

from flashgg.MicroAOD.flashggPhotons_cfi import flashggPhotons
from flashgg.MicroAOD.flashggRandomizedPhotonProducer_cff import flashggRandomizedPhotons
from flashgg.MicroAOD.flashggDiPhotons_cfi import flashggDiPhotons
from flashgg.MicroAOD.flashggMicroAODGenSequence_cff import *

RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService")
RandomNumberGeneratorService.flashggRandomizedPhotons = cms.PSet(
                          initialSeed = cms.untracked.uint32(16253245)
                          )

flashggDiphotonVertexSequence = cms.Sequence( flashggVertexMapUnique + flashggVertexMapNonUnique
                                             +flashggMicroAODGenSequence
                                             +flashggPhotons * flashggRandomizedPhotons * flashggDiPhotons
                                            )
