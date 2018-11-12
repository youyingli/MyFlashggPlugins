import os
import FWCore.ParameterSet.Config as cms

def includeflashggDiphoton(process):
    from flashgg.MicroAOD.flashggTkVtxMap_cfi import flashggVertexMapUnique,flashggVertexMapNonUnique
    setattr(process, 'flashggVertexMapUnique', flashggVertexMapUnique)
    setattr(process, 'flashggVertexMapNonUnique', flashggVertexMapNonUnique)

    process.load("flashgg.MicroAOD.flashggPhotons_cfi")
    process.load("flashgg.MicroAOD.flashggRandomizedPhotonProducer_cff")
    process.load("flashgg.MicroAOD.flashggDiPhotons_cfi")

def includeflashggLepton(process):
    process.load("flashgg.MicroAOD.flashggElectrons_cfi")
    process.load("flashgg.MicroAOD.flashggMuons_cfi")
    process.load("flashgg.MicroAOD.flashggLeptonSelectors_cff")

def includeflashggJet(process, isMC):
    from flashgg.MicroAOD.flashggTkVtxMap_cfi import flashggVertexMapForCHS
    from MyFlashggPlugins.flashggAnalysisNtuplizer.flashggJets_cfi import addFlashggPFCHSJets, maxJetCollections

    JetCollectionVInputTag = cms.VInputTag()
    for ivtx in range(0,maxJetCollections):
        addFlashggPFCHSJets (process       = process,
                             isData        = (not isMC),
                             vertexIndex   = ivtx,
                             useZeroVertex = False,
                             label         = '' + str(ivtx))
        JetCollectionVInputTag.append(cms.InputTag('flashggSelectedPFCHSJets' + str(ivtx)))
    flashggFinalJets = cms.EDProducer("FlashggVectorVectorJetCollector",
                                      inputTagJets = JetCollectionVInputTag
    )
    setattr(process, 'flashggVertexMapForCHS', flashggVertexMapForCHS)
    setattr(process, 'flashggFinalJets', flashggFinalJets)

def includePFMET(process, isMC):

    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
    runMetCorAndUncFromMiniAOD(process,
                         isData = (not isMC),
                         fixEE2017 = (not isMC),
                         fixEE2017Params = {'userawPt': True, 'PtThreshold':50.0, 'MinEtaThreshold':2.65, 'MaxEtaThreshold': 3.139},
                         # will produce new MET collection: slimmedMETsModifiedMET
                         postfix = "ModifiedMET",
                         )
    process.flashggMetsCorr = cms.EDProducer('FlashggMetProducer',
                                          verbose = cms.untracked.bool(False),
                                          metTag = cms.InputTag('slimmedMETsModifiedMET'),
                                          )

def includeSummer16EleID(process):
    from PhysicsTools.SelectorUtils.tools.vid_id_tools import DataFormat,switchOnVIDElectronIdProducer,setupAllVIDIdsInModule,setupVIDElectronSelection
    dataFormat = DataFormat.MiniAOD
    switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
    my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',
                     'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
                     'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff']
    for idmod in my_id_modules:
        setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
    process.flashggElectrons.effAreasConfigFile = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt")
    process.flashggElectrons.eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto")
    process.flashggElectrons.eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose")
    process.flashggElectrons.eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium")
    process.flashggElectrons.eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight")
    process.flashggElectrons.eleMVAMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90")
    process.flashggElectrons.eleMVATightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80")
    process.flashggElectrons.mvaValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values")

def includeFall17EleID(process):
    from PhysicsTools.SelectorUtils.tools.vid_id_tools import DataFormat,switchOnVIDElectronIdProducer,setupAllVIDIdsInModule,setupVIDElectronSelection
    dataFormat = DataFormat.MiniAOD
    switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
    my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff',
                     'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff',
                     'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V1_cff',
                     'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff']
    for idmod in my_id_modules:
        setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
    process.flashggElectrons.effAreasConfigFile = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_92X.txt")
    process.flashggElectrons.eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-veto")
    process.flashggElectrons.eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-loose")
    process.flashggElectrons.eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-medium")
    process.flashggElectrons.eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-tight")
    process.flashggElectrons.eleMVALooseIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wpLoose")
    process.flashggElectrons.eleMVAMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wp90")
    process.flashggElectrons.eleMVATightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wp80")
    process.flashggElectrons.eleMVALooseNoIsoIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wpLoose")
    process.flashggElectrons.eleMVAMediumNoIsoIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp90")
    process.flashggElectrons.eleMVATightNoIsoIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp80")
    process.flashggElectrons.mvaValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV1Values")
    process.flashggElectrons.mvaNoIsoValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Values")

def includeSummer16EGMPhoID(process):
    from PhysicsTools.SelectorUtils.tools.vid_id_tools import DataFormat,switchOnVIDPhotonIdProducer,setupAllVIDIdsInModule,setupVIDPhotonSelection
    dataFormat = DataFormat.MiniAOD
    switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD)
    my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff']
    for idmod in my_id_modules:
        setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

def includeFall17EGMPhoID(process):
    from PhysicsTools.SelectorUtils.tools.vid_id_tools import DataFormat,switchOnVIDPhotonIdProducer,setupAllVIDIdsInModule,setupVIDPhotonSelection
    dataFormat = DataFormat.MiniAOD
    switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD)
    my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1_cff']
    for idmod in my_id_modules:
        setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)
    process.flashggPhotons.effAreasConfigFile = cms.FileInPath("RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_TrueVtx.txt")
    process.flashggPhotons.egmMvaValuesMap = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRunIIFall17v1Values")

def includeflashggGenInfo(process):
    process.load("flashgg.MicroAOD.flashggMicroAODGenSequence_cff")

def includeflashggPDFs(process):
    process.load("flashgg/MicroAOD/flashggPDFWeightObject_cfi")

def includeHTXS(process):

    process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
    process.rivetProducerHTXS = cms.EDProducer('HTXSRivetProducer',
                                               HepMCCollection = cms.InputTag('myGenerator','unsmeared'),
                                               LHERunInfo = cms.InputTag('externalLHEProducer'),
                                               ProductionMode = cms.string('AUTO'),
                                               )
    process.mergedGenParticles = cms.EDProducer("MergedGenParticleProducer",
                                                inputPruned = cms.InputTag("prunedGenParticles"),
                                                inputPacked = cms.InputTag("packedGenParticles"),
                                                )
    process.myGenerator = cms.EDProducer("GenParticles2HepMCConverter",
                                         genParticles = cms.InputTag("mergedGenParticles"),
                                         genEventInfo = cms.InputTag("generator"),
                                         signalParticlePdgIds = cms.vint32(25), ## for the Higgs analysis
                                         )

# signal specific setting
def prepareSignal(process, filename):

    includeflashggDiphoton(process)
    includeflashggLepton(process)
    includeflashggJet(process, isMC = True)
    includePFMET(process, isMC = True)
    includeFall17EleID(process)
    includeFall17EGMPhoID(process)
    includeflashggGenInfo(process)
    includeHTXS(process)
    includeflashggPDFs(process)

    if filename.find("GluGlu") != -1:
        process.rivetProducerHTXS.ProductionMode = "GGF"
    if filename.find("THQ") != -1 or filename.find("THW") != -1:
        process.flashggPDFWeightObject.isStandardSample = False
        process.flashggPDFWeightObject.isThqSample = True

    process.flashggGenPhotonsExtra.defaultType = 1

# background specific setting
def prepareBackground(process, filename):

    includeflashggDiphoton(process)
    includeflashggLepton(process)
    includeflashggJet(process, isMC = True)
    includePFMET(process, isMC = True)
    includeFall17EleID(process)
    includeFall17EGMPhoID(process)
    includeflashggGenInfo(process)

    if filename.find("Sherpa") != -1:
        process.flashggGenPhotonsExtra.defaultType = 1

# data specific setting
def prepareData(process):

    includeflashggDiphoton(process)
    includeflashggLepton(process)
    includeflashggJet(process, isMC = False)
    includePFMET(process, isMC = False)
    includeFall17EleID(process)
    includeFall17EGMPhoID(process)

    from MyFlashggPlugins.flashggAnalysisNtuplizer.flashggJets2_cfi import maxJetCollections
    for vtx in range(0, maxJetCollections):
        delattr(process, "patJetGenJetMatchAK4PFCHSLeg%i"%vtx)
        delattr(process, "patJetFlavourAssociationAK4PFCHSLeg%i"%vtx)
        delattr(process, "patJetPartons%i"%vtx)
        delattr(process, "patJetPartonMatchAK4PFCHSLeg%i"%vtx)

def prepareflashggMicroAODTask(process, processType, filename):

    if processType == 'sig':
        prepareSignal(process, filename)
    elif processType == 'bkg':
        prepareBackground(process, filename)
    elif processType == 'data':
        prepareData(process, filename)
    else:
        raise Exception, "Please specify 'sig', 'bkg', 'data'"

    setattr( process, 'MicroAODTask', cms.Task() )
    getattr( process, 'MicroAODTask', cms.Task() ).add(*[getattr(process,prod) for prod in process.producers_()])
    getattr( process, 'MicroAODTask', cms.Task() ).add(*[getattr(process,filt) for filt in process.filters_()])
    return process.MicroAODTask

#################  H T X S  #########################
#process.rivetProducerHTXS.ProductionMode = "GGF", "VBF", "VH", "TTH", "TH"
