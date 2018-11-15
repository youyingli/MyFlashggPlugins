import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as opts

options = opts.VarParsing ('analysis')

options.register('Era',
                 "2017",
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'Era'
                 )

options.register('runMiniAOD',
                 False,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'runMiniAOD'
                 )

options.register('processType',
                 "sig",
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'processType'
                 )

options.register('filename',
                 "ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8",
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'filename'
                 )

options.register('doHTXS',
                 True,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'doHTXS'
                 )

options.register('doSystematics',
                 True,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'doSystematics'
                 )

options.parseArguments()

process = cms.Process("flashggAnalysisNtuplizerStd")

# geometry and global tag:
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
if options.processType == 'data':
    process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_ReReco_EOY17_v6')
else:
    process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v14')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(400) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1 )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.source = cms.Source ("PoolSource",
        fileNames = cms.untracked.vstring(
#'file:myMicroAODOutputFile_1.root'
#'file:myMicroAODOutputFile_4.root'
#'file:data.root'
#'/store/mc/RunIIFall17MiniAODv2/GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/90000/044F95FB-A342-E811-907F-5065F3816251.root'
'file:tth.root',
#'/store/group/phys_higgs/cmshgg/spigazzi/flashgg/RunIIFall17-3_2_0/RunIIFall17-3_2_0/DoubleEG/RunIIFall17-3_2_0-RunIIFall17-3_2_0-v0-Run2017F-09May2018-v1/181008_110542/0000/myMicroAODOutputFile_572.root'
#'/store/data/Run2017D/DoubleEG/MINIAOD/31Mar2018-v1/00000/002F7CD1-9D37-E811-A03E-B499BAABCF1A.root'
#'/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIIFall17-3_1_0/3_1_0/GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8_PSWeights/RunIIFall17-3_1_0-3_1_0-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/180605_202241/0000/myMicroAODOutputFile_3.root'
#'/store/mc/RunIIFall17MiniAODv2/GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/90000/EEC90FD7-A242-E811-AA78-EC0D9A8225FE.root'
#'/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIIFall17-3_1_0/3_1_0/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17-3_1_0-3_1_0-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180609_083905/0002/myMicroAODOutputFile_2040.root'
#'/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIIFall17-3_1_0/3_1_0/DoubleEG/RunIIFall17-3_1_0-3_1_0-v0-Run2017C-31Mar2018-v1/180606_155753/0000/myMicroAODOutputFile_217.root'
        )
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("DiPhotonNtuples.root")
)

#Modules builder
#**************************************************************
process.stdDiPhotonJetsSeq = cms.Sequence()

#---------------------------------------------------------------------------------------------
# Diphoton Trigger setting
# Data : Directly filter during processing 
# MC   : Store in Ntuple
#---------------------------------------------------------------------------------------------
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*",))
if options.processType == 'data':
    process.stdDiPhotonJetsSeq += process.hltHighLevel

#---------------------------------------------------------------------------------------------
# Dump ntuples from MiniAOD not microAOD. Very Slow!  Default is False
#---------------------------------------------------------------------------------------------
if options.runMiniAOD:
    process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService")
    process.RandomNumberGeneratorService.flashggRandomizedPhotons = cms.PSet(
              initialSeed = cms.untracked.uint32(16253245)
            )

    from MyFlashggPlugins.flashggAnalysisNtuplizer.prepareflashggMicroAODTask import prepareflashggMicroAODTask
    MicroAODTask = prepareflashggMicroAODTask(process, options.processType, options.filename)
    process.stdDiPhotonJetsSeq.associate( MicroAODTask )

#---------------------------------------------------------------------------------------------
# DiPhoton preselection setting
#---------------------------------------------------------------------------------------------
process.load("flashgg/Taggers/flashggUpdatedIdMVADiPhotons_cfi")
process.load("flashgg/Taggers/flashggPreselectedDiPhotons_cfi")

#---------------------------------------------------------------------------------------------
# Jets unpack setting
#---------------------------------------------------------------------------------------------
from MyFlashggPlugins.flashggAnalysisNtuplizer.flashggJets_cfi import maxJetCollections
flashggUnpackedJets = cms.EDProducer("FlashggVectorVectorJetUnpacker",
                                     JetsTag = cms.InputTag("flashggFinalJets"),
                                     NCollections = cms.uint32(maxJetCollections)
                                    )

setattr(process, 'flashggUnpackedJets', flashggUnpackedJets)

UnpackedJetCollectionVInputTag = cms.VInputTag()
for i in range(0,maxJetCollections):
    UnpackedJetCollectionVInputTag.append(cms.InputTag('flashggUnpackedJets',str(i)))

#---------------------------------------------------------------------------------------------
# Merge DiPhoton and Jets modules
#---------------------------------------------------------------------------------------------
process.basicSeq = cms.Sequence(process.flashggUpdatedIdMVADiPhotons
                               *process.flashggPreselectedDiPhotons
                               *process.flashggUnpackedJets
                                )

process.stdDiPhotonJetsSeq += process.basicSeq

#---------------------------------------------------------------------------------------------
# Systematics setting (To Do)
#---------------------------------------------------------------------------------------------
from MyFlashggPlugins.flashggAnalysisNtuplizer.prepareflashggDiPhotonSystematicsTask import prepareflashggDiPhotonSystematicsTask, getDiPhotonSystematicsList
diphotonSystematicsTask = prepareflashggDiPhotonSystematicsTask(process, options.processType, True)
diphosystname = ['']
diphoton = [cms.InputTag('flashggPreselectedDiPhotons')]

if options.doSystematics:
    for syst in getDiPhotonSystematicsList():
        diphosystname.append(syst)
        diphoton.append(cms.InputTag('flashggPreselectedDiPhotons' + syst))

#---------------------------------------------------------------------------------------------
# Customize your own settings based on exist modules
#---------------------------------------------------------------------------------------------
from MyFlashggPlugins.flashggAnalysisNtuplizer.myflashggCustomize import myflashggCustomize
myflashggCustomize(process, options.runMiniAOD)

#---------------------------------------------------------------------------------------------
# Main Ntuplizer
#---------------------------------------------------------------------------------------------
from flashgg.Taggers.flashggTags_cff import HTXSInputTags
process.flashggNtuples = cms.EDAnalyzer('flashggAnaTreeMerge',
                                       diphosystnames  = cms.vstring(diphosystname),
                                       diphotons       = cms.VInputTag(diphoton),
                                       nondiphosetting = cms.PSet(
                                           inputTagJets            = UnpackedJetCollectionVInputTag,
                                           ElectronTag             = cms.InputTag('flashggSelectedElectrons'),
                                           MuonTag                 = cms.InputTag('flashggSelectedMuons'),
                                           MetTag                  = cms.InputTag('flashggMetsCorr'),
                                           VertexTag               = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                           BeamSpotTag             = cms.InputTag('offlineBeamSpot'),
                                           RhoTag                  = cms.InputTag('fixedGridRhoFastjetAll'),
                                           GenParticleTag          = cms.InputTag('flashggPrunedGenParticles'),
                                           GenEventInfo            = cms.InputTag('generator'),
                                           PileUpTag               = cms.InputTag('slimmedAddPileupInfo'),
                                           TriggerTag              = cms.InputTag('TriggerResults::HLT'),
                                           MetTriggerTag           = cms.InputTag('TriggerResults::PAT'),
                                           pathName                = cms.string("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90"),
                                           isMiniAOD               = cms.bool(options.runMiniAOD),
                                           doHTXS                  = cms.bool(options.doHTXS),
                                           HTXSTags                = HTXSInputTags
                                       )
)

process.stdDiPhotonJetsSeq += process.flashggNtuples

#---------------------------------------------------------------------------------------------
# Final Path to run
#---------------------------------------------------------------------------------------------
process.p = cms.Path( process.stdDiPhotonJetsSeq, diphotonSystematicsTask)
