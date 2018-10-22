import FWCore.ParameterSet.Config as cms
from SimGeneral.HepPDTESSource.pythiapdt_cfi import HepPDTESSource

rivetProducerHTXS = cms.EDProducer('HTXSRivetProducer',
                                   HepMCCollection = cms.InputTag('myGenerator','unsmeared'),
                                   LHERunInfo = cms.InputTag('externalLHEProducer'),
                                   ProductionMode = cms.string('AUTO'),
                                   )

mergedGenParticles = cms.EDProducer("MergedGenParticleProducer",
                                    inputPruned = cms.InputTag("prunedGenParticles"),
                                    inputPacked = cms.InputTag("packedGenParticles"),
                                    )

myGenerator = cms.EDProducer("GenParticles2HepMCConverter",
                             genParticles = cms.InputTag("mergedGenParticles"),
                             genEventInfo = cms.InputTag("generator"),
                             signalParticlePdgIds = cms.vint32(25), ## for the Higgs analysis
                             )
