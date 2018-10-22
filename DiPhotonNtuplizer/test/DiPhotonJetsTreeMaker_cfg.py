import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as opts

options = opts.VarParsing ('analysis')
options.register('runMiniAOD',
                 False,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'runMiniAOD'
                 )

options.register('isData',
                 False,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'isData'
                 )

options.register('doHTXS',
                 False,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'doHTXS'
                 )

options.parseArguments()

process = cms.Process("DiPhotonJetsTreeMaker")

# geometry and global tag:
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
if options.isData:
    process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_ReReco_EOY17_v6')
else:
    process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v14')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 100 )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.source = cms.Source ("PoolSource",
        fileNames = cms.untracked.vstring(
'file:myMicroAODOutputFile_1.root'
#'file:data.root'
#'/store/mc/RunIIFall17MiniAODv2/GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/90000/044F95FB-A342-E811-907F-5065F3816251.root'
#'/store/data/Run2017D/DoubleEG/MINIAOD/31Mar2018-v1/00000/002F7CD1-9D37-E811-A03E-B499BAABCF1A.root'
#'/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIIFall17-3_1_0/3_1_0/GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8_PSWeights/RunIIFall17-3_1_0-3_1_0-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/180605_202241/0000/myMicroAODOutputFile_3.root'
#'/store/mc/RunIIFall17MiniAODv2/GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/90000/EEC90FD7-A242-E811-AA78-EC0D9A8225FE.root'
#'/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIIFall17-3_1_0/3_1_0/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17-3_1_0-3_1_0-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180609_083905/0002/myMicroAODOutputFile_2040.root'
#'/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIIFall17-3_1_0/3_1_0/DoubleEG/RunIIFall17-3_1_0-3_1_0-v0-Run2017C-31Mar2018-v1/180606_155753/0000/myMicroAODOutputFile_217.root'
        )
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("DiPhotonJetsTree.root")
)

#Sequence builder
#**************************************************************
process.stdDiPhotonJetsSeq = cms.Sequence()
#-----------------------------------------------------------------------------------------------
from myflashggPlugins.vbf_analysis.flashggJets_cfi import maxJetCollections

if options.runMiniAOD:
    from myflashggPlugins.vbf_analysis.prepareflashggMicroAODSequence_cff import prepareMicroAODTask
    MicroAODTask = prepareMicroAODTask(process, options.isData)
    process.stdDiPhotonJetsSeq.associate( MicroAODTask )

flashggUnpackedJets = cms.EDProducer("FlashggVectorVectorJetUnpacker",
                                     JetsTag = cms.InputTag("flashggFinalJets"),
                                     NCollections = cms.uint32(maxJetCollections)
                                    )

setattr(process, 'flashggUnpackedJets', flashggUnpackedJets)

UnpackedJetCollectionVInputTag = cms.VInputTag()
for i in range(0,maxJetCollections):
    UnpackedJetCollectionVInputTag.append(cms.InputTag('flashggUnpackedJets',str(i)))

#process.load("flashggDiphotonVtxPlugins/MVATraining/flashggDiphotonVertexSequence_cff")
#process.flashggDiPhotons.vertexIdMVAweightfile = cms.FileInPath("flashggDiphotonVtxPlugins/Validation/data/TMVAClassification_BDTVtxId_SL_2017V3.xml")
#process.flashggDiPhotons.vertexProbMVAweightfile = cms.FileInPath("flashggDiphotonVtxPlugins/Validation/data/TMVAClassification_BDTVtxProb_SL_2017V3.xml")
#process.flashggDiPhotons.nVtxSaveInfo = cms.untracked.uint32(999)
#process.flashggDiPhotons.useSingleLeg = cms.bool(True)
#process.flashggDiPhotons.useZerothVertexFromMicro = cms.bool(False)
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

process.load("flashgg/Taggers/flashggPreselectedDiPhotons_cfi")
process.load("flashgg/Taggers/flashggUpdatedIdMVADiPhotons_cfi")

process.basicSeq = cms.Sequence(process.flashggUpdatedIdMVADiPhotons
                               *process.flashggPreselectedDiPhotons
                               *process.flashggUnpackedJets
                                )

process.stdDiPhotonJetsSeq += process.basicSeq

from flashgg.Taggers.flashggTags_cff import HTXSInputTags
process.commissioning = cms.EDAnalyzer('DiPhotonJetsTreeMaker',
                                       DiPhotonTag             = cms.InputTag('flashggPreselectedDiPhotons'),
                                       inputTagJets            = UnpackedJetCollectionVInputTag,
                                       VertexTag               = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                       VertexCandidateMapTagDz = cms.InputTag('flashggVertexMapUnique'),
                                       BeamSpotTag             = cms.InputTag('offlineBeamSpot'),
                                       RhoTag                  = cms.InputTag('fixedGridRhoFastjetAll'),
                                       #GenParticleTag          = cms.InputTag('prunedGenParticles'),
                                       GenParticleTag          = cms.InputTag('flashggPrunedGenParticles'),
                                       GenEventInfo            = cms.InputTag('generator'),
                                       PileUpTag               = cms.InputTag('slimmedAddPileupInfo'),
                                       TriggerTag              = cms.InputTag('TriggerResults::HLT'),
                                       pathName                = cms.string("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90"),
                                       doHTXS                  = cms.bool(options.doHTXS),
                                       HTXSTags                = HTXSInputTags
)

process.stdDiPhotonJetsSeq += process.commissioning

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*",))
if options.isData:
    process.stdDiPhotonJetsSeq += process.hltHighLevel

process.p = cms.Path( process.stdDiPhotonJetsSeq )

#process.p = cms.Path(process.flashggUpdatedIdMVADiPhotons
#                    *process.flashggPreselectedDiPhotons
#                    *process.flashggUnpackedJets
#                    #*process.hltHighLevel
#                    *process.commissioning#,  MicroAODTask
#                    )
