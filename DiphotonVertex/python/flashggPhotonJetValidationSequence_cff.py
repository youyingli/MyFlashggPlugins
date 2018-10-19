import FWCore.ParameterSet.Config as cms
from flashggDiphotonVtxPlugins.Validation.flashggTkVtxMapValidation_cfi import flashggVertexMapUnique,flashggVertexMapNonUnique
from flashgg.MicroAOD.flashggPhotons_cfi import flashggPhotons
from flashggDiphotonVtxPlugins.Validation.flashggPhotonJet_cfi import flashggPhotonJet
from flashgg.MicroAOD.flashggMicroAODGenSequence_cff import *

flashggPhotonJetValidationSequence = cms.Sequence( flashggVertexMapUnique + flashggVertexMapNonUnique
                                                  +flashggMicroAODGenSequence
                                                  +flashggPhotons
                                                  +flashggPhotonJet
                                                 )
