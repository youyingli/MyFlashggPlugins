import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("PhotonJetValidationTreeMaker")

# geometry and global tag:
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_ReReco_EOY17_v2')
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v12')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 100 )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.source = cms.Source ("PoolSource",
        fileNames = cms.untracked.vstring(
'/store/mc/RunIIFall17MiniAOD/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v2/10000/28BD4EBE-67F8-E711-A645-5065F37D4131.root'
        )
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("PhotonJetTree.root")
)


#Sequence builder
#**************************************************************
process.load("flashggDiphotonVtxPlugins.Validation.flashggPhotonJetValidationSequence_cff")
#process.flashggPhotonJetValidationSequence.remove(process.flashggMicroAODGenSequence)

process.flashggPhotonJet.vertexIdMVAweightfile = cms.FileInPath("flashggDiphotonVtxPlugins/Validation/data/TMVAClassification_BDTVtxId_SL_2017.xml")
process.flashggPhotonJet.vertexProbMVAweightfile = cms.FileInPath("flashggDiphotonVtxPlugins/Validation/data/TMVAClassification_BDTVtxProb_SL_2017.xml")
process.flashggPhotonJet.useSingleLeg=cms.bool(True)
process.flashggPhotonJet.sigma1Pix               = cms.double( 0.00800379 )
process.flashggPhotonJet.sigma1Tib               = cms.double( 0.502127   )
process.flashggPhotonJet.sigma1Tob               = cms.double( 5.5138     )
process.flashggPhotonJet.sigma1PixFwd            = cms.double( 0.0318172  )
process.flashggPhotonJet.sigma1Tid               = cms.double( 0.325117   )
process.flashggPhotonJet.sigma1Tec               = cms.double( 1.19907    )
process.flashggPhotonJet.sigma2Pix               = cms.double( 0.0171381  )
process.flashggPhotonJet.sigma2Tib               = cms.double( 0.282616   )
process.flashggPhotonJet.sigma2Tob               = cms.double( 3.5737     )
process.flashggPhotonJet.sigma2PixFwd            = cms.double( 0.0923745  )
process.flashggPhotonJet.sigma2Tid               = cms.double( 0.355705   )
process.flashggPhotonJet.sigma2Tec               = cms.double( 0.863342   )
process.flashggPhotonJet.singlelegsigma1Pix      = cms.double( 0.00879849 )
process.flashggPhotonJet.singlelegsigma1Tib      = cms.double( 1.37155    )
process.flashggPhotonJet.singlelegsigma1Tob      = cms.double( 2.7242     )
process.flashggPhotonJet.singlelegsigma1PixFwd   = cms.double( 0.0596455  )
process.flashggPhotonJet.singlelegsigma1Tid      = cms.double( 0.479279   )
process.flashggPhotonJet.singlelegsigma1Tec      = cms.double( 2.02211    )
process.flashggPhotonJet.singlelegsigma2Pix      = cms.double( 0.0224474  )
process.flashggPhotonJet.singlelegsigma2Tib      = cms.double( 0.594662   )
process.flashggPhotonJet.singlelegsigma2Tob      = cms.double( 0.433137   )
process.flashggPhotonJet.singlelegsigma2PixFwd   = cms.double( 0.137922   )
process.flashggPhotonJet.singlelegsigma2Tid      = cms.double( 0.421378   )
process.flashggPhotonJet.singlelegsigma2Tec      = cms.double( 0.977421   )


process.commissioning = cms.EDAnalyzer('PhotonJetValidationTreeMaker',
                                       PhotonJetTag=cms.InputTag('flashggPhotonJet'),
                                       vertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                       GenEventInfo=cms.InputTag('generator'),
                                       BeamSpotTag=cms.InputTag('offlineBeamSpot'),
                                       PileUpTag=cms.InputTag('slimmedAddPileupInfo'),
#                                       evWeight = cms.double(1.00000) #Data
                                       evWeight = cms.double(1.00000) #GJets_HT-40To100
#                                       evWeight = cms.double(0.44768) #GJets_HT-100To200
#                                       evWeight = cms.double(0.11234) #GJets_HT-200To400
#                                       evWeight = cms.double(0.01332) #GJets_HT-400To600
#                                       evWeight = cms.double(0.00451) #GJets_HT-600ToInf
)

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Photon50_v*"))

from PhysicsTools.SelectorUtils.tools.vid_id_tools import DataFormat,switchOnVIDPhotonIdProducer,setupAllVIDIdsInModule,setupVIDPhotonSelection
dataFormat = DataFormat.MiniAOD
switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD)
#my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff']
my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1_cff']
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)
process.flashggPhotons.effAreasConfigFile = cms.FileInPath("RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_TrueVtx.txt")
process.flashggPhotons.egmMvaValuesMap = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRunIIFall17v1Values")

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService")
process.RandomNumberGeneratorService.flashggRandomizedPhotons = cms.PSet(
                  initialSeed = cms.untracked.uint32(16253245)
                          )

process.p = cms.Path(process.flashggPhotonJetValidationSequence
                    *process.egmPhotonIDSequence
                    *process.hltHighLevel
                    *process.commissioning
                    )
