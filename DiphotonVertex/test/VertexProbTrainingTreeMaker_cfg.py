import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as opts

options = opts.VarParsing ('analysis')
options.register('isZeroVertex',
                 False,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'isZeroVertex'
                 )

options.parseArguments()

process = cms.Process("VertexProbTrainingTreeMaker")

# geometry and global tag:
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v18')
#process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v17')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 100 )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.source = cms.Source ("PoolSource",
        fileNames = cms.untracked.vstring(
            options.inputFiles
#'/store/mc/RunIIFall17MiniAODv2/GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/90000/044F95FB-A342-E811-907F-5065F3816251.root'
#'/store/mc/RunIIAutumn18MiniAOD/GluGluHToGG_M-125_13TeV_powheg_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/FF3573F4-A9DA-7B4C-A1BE-E893C22B75E3.root'
#'/store/mc/RunIIFall17MiniAODv2/ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/307FE799-6B42-E811-8A34-0025905A48D6.root'
        )
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.outputFile
                                       #"VertexProbTrainingTree.root"
                                       )
)

#Sequence builder
#**************************************************************
process.load("MyFlashggPlugins.DiphotonVertex.flashggDiphotonVertexSequence_cff")
process.flashggDiPhotons.vertexIdMVAweightfile = cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_BDTVtxId_SL_2016.xml")
process.flashggDiPhotons.vertexProbMVAweightfile = cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_BDTVtxProb_SL_2016.xml")
process.flashggDiPhotons.nVtxSaveInfo = cms.untracked.uint32(999)
process.flashggDiPhotons.useSingleLeg = cms.bool(True)
process.flashggDiPhotons.useZerothVertexFromMicro = cms.bool(options.isZeroVertex)

process.commissioning = cms.EDAnalyzer('VertexProbTrainingTreeMaker',
                                       DiPhotonTag             = cms.InputTag('flashggPreselectedDiPhotons'),
                                       VertexTag               = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                       VertexCandidateMapTagDz = cms.InputTag('flashggVertexMapUnique'),
                                       BeamSpotTag             = cms.InputTag('offlineBeamSpot'),
                                       rhoTag                  = cms.InputTag('fixedGridRhoFastjetAll'),
                                       GenParticleTag          = cms.InputTag('prunedGenParticles'),
                                       GenEventInfoTag         = cms.InputTag('generator'),
                                       PileUpTag               = cms.InputTag('slimmedAddPileupInfo'),
                                       isZeroVertex            = cms.bool(options.isZeroVertex)
)

from PhysicsTools.SelectorUtils.tools.vid_id_tools import DataFormat,switchOnVIDPhotonIdProducer,setupAllVIDIdsInModule,setupVIDPhotonSelection
dataFormat = DataFormat.MiniAOD
switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD)
my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V2_cff']
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)
process.flashggPhotons.photonIdMVAweightfile_EB = cms.FileInPath("flashgg/MicroAOD/data/HggPhoId_94X_barrel_BDT_woisocorr.weights.xml")
process.flashggPhotons.photonIdMVAweightfile_EE = cms.FileInPath("flashgg/MicroAOD/data/HggPhoId_94X_endcap_BDT_woisocorr.weights.xml")
process.flashggPhotons.effAreasConfigFile = cms.FileInPath("RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_TrueVtx.txt")
process.flashggPhotons.egmMvaValuesMap = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRunIIFall17v2Values")
process.flashggPhotons.is2017 = cms.bool(True)

process.load("flashgg/Taggers/flashggTagSequence_cfi")

process.p = cms.Path(process.flashggDiphotonVertexSequence
                    *process.egmPhotonIDSequence
                    *process.flashggUpdatedIdMVADiPhotons
                    *process.flashggPreselectedDiPhotons
                    *process.commissioning
                    )
