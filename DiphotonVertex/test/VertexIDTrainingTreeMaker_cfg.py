import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("VertexIDTrainingTreeMaker")

# geometry and global tag:
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v14')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 100 )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.source = cms.Source ("PoolSource",
        fileNames = cms.untracked.vstring(
'/store/mc/RunIIFall17MiniAODv2/GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/90000/044F95FB-A342-E811-907F-5065F3816251.root'
        )
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("VertexIDTrainingTree.root")
)

#Sequence builder
#**************************************************************
process.load("flashggDiphotonVtxPlugins/MVATraining/flashggDiphotonVertexSequence_cff")
process.flashggDiPhotons.nVtxSaveInfo = cms.untracked.uint32(999)
process.flashggDiPhotons.useSingleLeg = cms.bool(True)

#Cut Based Photon ID
#process.flashggDiPhotons.sigma1Pix               = cms.double( 0.00800379 )
#process.flashggDiPhotons.sigma1Tib               = cms.double( 0.502127   )
#process.flashggDiPhotons.sigma1Tob               = cms.double( 5.5138     )
#process.flashggDiPhotons.sigma1PixFwd            = cms.double( 0.0318172  )
#process.flashggDiPhotons.sigma1Tid               = cms.double( 0.325117   )
#process.flashggDiPhotons.sigma1Tec               = cms.double( 1.19907    )
#process.flashggDiPhotons.sigma2Pix               = cms.double( 0.0171381  )
#process.flashggDiPhotons.sigma2Tib               = cms.double( 0.282616   )
#process.flashggDiPhotons.sigma2Tob               = cms.double( 3.5737     )
#process.flashggDiPhotons.sigma2PixFwd            = cms.double( 0.0923745  )
#process.flashggDiPhotons.sigma2Tid               = cms.double( 0.355705   )
#process.flashggDiPhotons.sigma2Tec               = cms.double( 0.863342   )
#process.flashggDiPhotons.singlelegsigma1Pix      = cms.double( 0.00879849 )
#process.flashggDiPhotons.singlelegsigma1Tib      = cms.double( 1.37155    )
#process.flashggDiPhotons.singlelegsigma1Tob      = cms.double( 2.7242     )
#process.flashggDiPhotons.singlelegsigma1PixFwd   = cms.double( 0.0596455  )
#process.flashggDiPhotons.singlelegsigma1Tid      = cms.double( 0.479279   )
#process.flashggDiPhotons.singlelegsigma1Tec      = cms.double( 2.02211    )
#process.flashggDiPhotons.singlelegsigma2Pix      = cms.double( 0.0224474  )
#process.flashggDiPhotons.singlelegsigma2Tib      = cms.double( 0.594662   )
#process.flashggDiPhotons.singlelegsigma2Tob      = cms.double( 0.433137   )
#process.flashggDiPhotons.singlelegsigma2PixFwd   = cms.double( 0.137922   )
#process.flashggDiPhotons.singlelegsigma2Tid      = cms.double( 0.421378   )
#process.flashggDiPhotons.singlelegsigma2Tec      = cms.double( 0.977421   )

#Diphoton preselection
#process.flashggDiPhotons.sigma1Pix               = cms.double( 0.00802885 )
#process.flashggDiPhotons.sigma1Tib               = cms.double( 0.518602   )
#process.flashggDiPhotons.sigma1Tob               = cms.double( 5.54238    )
#process.flashggDiPhotons.sigma1PixFwd            = cms.double( 0.0343765  )
#process.flashggDiPhotons.sigma1Tid               = cms.double( 0.327043   )
#process.flashggDiPhotons.sigma1Tec               = cms.double( 1.22594    )
#process.flashggDiPhotons.sigma2Pix               = cms.double( 0.0159439  )
#process.flashggDiPhotons.sigma2Tib               = cms.double( 0.285931   )
#process.flashggDiPhotons.sigma2Tob               = cms.double( 3.74779    )
#process.flashggDiPhotons.sigma2PixFwd            = cms.double( 0.0914724  )
#process.flashggDiPhotons.sigma2Tid               = cms.double( 0.342842   )
#process.flashggDiPhotons.sigma2Tec               = cms.double( 0.875536   )
#process.flashggDiPhotons.singlelegsigma1Pix      = cms.double( 0.00886973 )
#process.flashggDiPhotons.singlelegsigma1Tib      = cms.double( 1.36891    )
#process.flashggDiPhotons.singlelegsigma1Tob      = cms.double( 2.69513    )
#process.flashggDiPhotons.singlelegsigma1PixFwd   = cms.double( 0.064971   )
#process.flashggDiPhotons.singlelegsigma1Tid      = cms.double( 0.457494   )
#process.flashggDiPhotons.singlelegsigma1Tec      = cms.double( 1.98512    )
#process.flashggDiPhotons.singlelegsigma2Pix      = cms.double( 0.0218862  )
#process.flashggDiPhotons.singlelegsigma2Tib      = cms.double( 0.581801   )
#process.flashggDiPhotons.singlelegsigma2Tob      = cms.double( 0.42568    )
#process.flashggDiPhotons.singlelegsigma2PixFwd   = cms.double( 0.139077   )
#process.flashggDiPhotons.singlelegsigma2Tid      = cms.double( 0.395489   )
#process.flashggDiPhotons.singlelegsigma2Tec      = cms.double( 0.958625   )

process.flashggDiPhotons.sigma1Pix               = cms.double( 0.00804137 )
process.flashggDiPhotons.sigma1Tib               = cms.double( 0.507857   )
process.flashggDiPhotons.sigma1Tob               = cms.double( 5.63449    )
process.flashggDiPhotons.sigma1PixFwd            = cms.double( 0.034551   )
process.flashggDiPhotons.sigma1Tid               = cms.double( 0.326273   )
process.flashggDiPhotons.sigma1Tec               = cms.double( 1.19907    )
process.flashggDiPhotons.sigma2Pix               = cms.double( 0.0160633  )
process.flashggDiPhotons.sigma2Tib               = cms.double( 0.282202   )
process.flashggDiPhotons.sigma2Tob               = cms.double( 3.95819    )
process.flashggDiPhotons.sigma2PixFwd            = cms.double( 0.0920137  )
process.flashggDiPhotons.sigma2Tid               = cms.double( 0.342842   )
process.flashggDiPhotons.sigma2Tec               = cms.double( 0.875536   )
process.flashggDiPhotons.singlelegsigma1Pix      = cms.double( 0.00956435 )
process.flashggDiPhotons.singlelegsigma1Tib      = cms.double( 1.38078    )
process.flashggDiPhotons.singlelegsigma1Tob      = cms.double( 2.60343    )
process.flashggDiPhotons.singlelegsigma1PixFwd   = cms.double( 0.0835342  )
process.flashggDiPhotons.singlelegsigma1Tid      = cms.double( 0.471549   )
process.flashggDiPhotons.singlelegsigma1Tec      = cms.double( 1.99499    )
process.flashggDiPhotons.singlelegsigma2Pix      = cms.double( 0.0288075  )
process.flashggDiPhotons.singlelegsigma2Tib      = cms.double( 0.58634    )
process.flashggDiPhotons.singlelegsigma2Tob      = cms.double( 0.400822   )
process.flashggDiPhotons.singlelegsigma2PixFwd   = cms.double( 0.171393   )
process.flashggDiPhotons.singlelegsigma2Tid      = cms.double( 0.403523   )
process.flashggDiPhotons.singlelegsigma2Tec      = cms.double( 0.958625   )

process.commissioning = cms.EDAnalyzer('VertexIDTrainingTreeMaker',
                                       #DiPhotonTag             = cms.InputTag('flashggDiPhotons'),
                                       DiPhotonTag             = cms.InputTag('flashggPreselectedDiPhotons'),
                                       VertexTag               = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                       VertexCandidateMapTagDz = cms.InputTag('flashggVertexMapUnique'),
                                       BeamSpotTag             = cms.InputTag('offlineBeamSpot'),
                                       rhoTag                  = cms.InputTag('fixedGridRhoFastjetAll'),
                                       GenParticleTag          = cms.InputTag('prunedGenParticles'),
                                       GenEventInfo            = cms.InputTag('generator')
)

from PhysicsTools.SelectorUtils.tools.vid_id_tools import DataFormat,switchOnVIDPhotonIdProducer,setupAllVIDIdsInModule,setupVIDPhotonSelection
dataFormat = DataFormat.MiniAOD
switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD)
my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1_cff']
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)
process.flashggPhotons.effAreasConfigFile = cms.FileInPath("RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_TrueVtx.txt")
process.flashggPhotons.egmMvaValuesMap = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRunIIFall17v1Values")

process.load("flashgg/Taggers/flashggTagSequence_cfi")

process.p = cms.Path(process.flashggDiphotonVertexSequence
                    *process.egmPhotonIDSequence
                    *process.flashggUpdatedIdMVADiPhotons
                    *process.flashggPreselectedDiPhotons
                    *process.commissioning
                    )
