# this is based on the jet 74X recipe for jet re-clustering
# twiki source : https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015#Advanced_topics_re_clustering_ev

import FWCore.ParameterSet.Config as cms

#from RecoJets.JetProducers.PileupJetIDParams_cfi import cutbased_new as pu_jetid
from PhysicsTools.PatAlgos.tools.jetTools        import addJetCollection
from CondCore.DBCommon.CondDBSetup_cfi import *

import os

flashggBTag = 'pfCombinedInclusiveSecondaryVertexV2BJetTags'

maxJetCollections = 12
#qgDatabaseVersion = 'v1' # check https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion
if os.environ["CMSSW_VERSION"].count("CMSSW_9"):
  qgDatabaseVersion = 'cmssw8020_v2'
else:
  qgDatabaseVersion = '80X'

def addFlashggPFCHSJets(process, 
                        isData,
                        vertexIndex = 0, 
                        #doQGTagging = True, 
                        label ='',
                        debug = False):

    from flashgg.MicroAOD.flashggTkVtxMap_cfi import flashggVertexMapForCHS
    setattr(process, 'flashggVertexMapForCHS', flashggVertexMapForCHS)

  
    setattr(process, 'selectedMuons' + label, cms.EDFilter("CandPtrSelector", 
                                                           src = cms.InputTag("slimmedMuons"), 
                                                           cut = cms.string('''abs(eta)<2.5 && pt>10. &&
                                                         (pfIsolationR04().sumChargedHadronPt+
                                                         max(0.,pfIsolationR04().sumNeutralHadronEt+
                                                         pfIsolationR04().sumPhotonEt-
                                                         0.50*pfIsolationR04().sumPUPt))/pt < 0.20 && 
                                                         (isPFMuon && (isGlobalMuon || isTrackerMuon) )''')))
    
    setattr(process, 'selectedElectrons' + label, cms.EDFilter("CandPtrSelector", 
                                                               src = cms.InputTag("slimmedElectrons"), 
                                                               cut = cms.string('''abs(eta)<2.5 && pt>20. &&
                                                             gsfTrack.isAvailable() &&
                                                             gsfTrack.hitPattern().numberOfLostHits(\'MISSING_INNER_HITS\') < 2 &&
                                                             (pfIsolationVariables().sumChargedHadronPt+
                                                             max(0.,pfIsolationVariables().sumNeutralHadronEt+
                                                             pfIsolationVariables().sumPhotonEt-
                                                             0.5*pfIsolationVariables().sumPUPt))/pt < 0.15''')))
    
    # Simple producer which just removes the Candidates which
    # don't come from the legacy vertex according to the Flashgg Vertex Map
    setattr(process,'flashggCHSLegacyVertexCandidates' + label,
            cms.EDProducer('FlashggMultiCHSLegacyVertexCandProducer',
                           PFCandidatesTag       = cms.InputTag('packedPFCandidates'),
                           DiPhotonTag           = cms.InputTag('flashggDiPhotons'),
                           VertexCandidateMapTag = cms.InputTag("flashggVertexMapForCHS"),
                           VertexTag             = cms.InputTag('offlineSlimmedPrimaryVertices'),
                           vertexIndex           = cms.uint32(vertexIndex),
                           debug                 = cms.untracked.bool(debug)
                           )
            )
    
    setattr(process, 'pfCHSLeg' + label, cms.EDFilter("CandPtrSelector", 
                                                      src = cms.InputTag('flashggCHSLegacyVertexCandidates' + label), 
                                                      cut = cms.string('')))
    
    # then remove the previously selected muons
    setattr(process, 'pfNoMuonCHSLeg' + label,  cms.EDProducer("CandPtrProjector", 
                                                               src  = cms.InputTag("pfCHSLeg" + label), 
                                                               veto = cms.InputTag("selectedMuons" + label)))
    # then remove the previously selected electrons
    setattr(process, 'pfNoElectronsCHSLeg' + label,  cms.EDProducer("CandPtrProjector", 
                                                                    src  = cms.InputTag("pfNoMuonCHSLeg" + label), 
                                                                    veto = cms.InputTag("selectedElectrons" + label)))
  
   
    #Import RECO jet producer for ak4 PF and GEN jet
    from RecoJets.JetProducers.ak4PFJets_cfi  import ak4PFJets
    setattr(process, 'ak4PFJetsCHSLeg' + label, ak4PFJets.clone ( src = 'pfNoElectronsCHSLeg' + label, doAreaFastjet = True))
  
    if isData:
        JECs = ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual']
    else:
        JECs = ['L1FastJet', 'L2Relative', 'L3Absolute']
  
    # NOTE: this is the 74X recipe for the jet clustering
    addJetCollection(
        process,
        postfix        = label,
        labelName      = 'AK4PFCHSLeg',
        jetSource      = cms.InputTag('ak4PFJetsCHSLeg' + label),
        pvSource       = cms.InputTag('offlineSlimmedPrimaryVertices'),
        pfCandidates   = cms.InputTag('packedPFCandidates'),
        svSource       = cms.InputTag('slimmedSecondaryVertices'),
        elSource       = cms.InputTag("slimmedElectrons"),
        muSource       = cms.InputTag("slimmedMuons"),
        runIVF         = True,
        btagDiscriminators = [ flashggBTag ],
        jetCorrections = ('AK4PFchs', JECs, 'None'),
        genJetCollection = cms.InputTag('slimmedGenJets'),
        genParticles     = cms.InputTag('prunedGenParticles'),
        # jet param
        algo = 'AK', rParam = 0.4
        #btagInfos =  ['pfImpactParameterTagInfos', 'pfSecondaryVertexTagInfos', 'pfInclusiveSecondaryVertexFinderTagInfos', 'softPFMuonsTagInfos', 'softPFElectronsTagInfos'] #Extra Btagging Info
    )
  
    #Recalculate btagging info
    getattr( process, 'patJetsAK4PFCHSLeg' + label).addTagInfos = True
    
    #adjust PV used for Jet Corrections
    getattr(process, 'patJetCorrFactorsAK4PFCHSLeg' + label).primaryVertices = "offlineSlimmedPrimaryVertices"
  
  
    if not hasattr(process,"QGPoolDBESSource"):
        process.QGPoolDBESSource = cms.ESSource("PoolDBESSource",
                                                CondDBSetup,
                                                toGet = cms.VPSet(
                                                        cms.PSet(
                                                            record = cms.string('QGLikelihoodRcd'),
                                                            tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_AK4PFchs'),
                                                            label  = cms.untracked.string('QGL_AK4PFchs')
                                                        )
                                                ),
                                                connect = cms.string('sqlite_fip:flashggDiphotonVtxPlugins/MVATraining/data/QGL_'+qgDatabaseVersion+'.db') 
        )
        process.es_prefer_qg = cms.ESPrefer('PoolDBESSource','QGPoolDBESSource')
    
    for type in ['AK4PFchs']:#,'AK4PFchs_antib']:
        process.QGPoolDBESSource.toGet.extend(cms.VPSet(cms.PSet(
              record = cms.string('QGLikelihoodRcd'),
              tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_'+type),
              label  = cms.untracked.string('QGL_'+type)
              )))
    
    from RecoJets.JetProducers.QGTagger_cfi import QGTagger
    setattr( process, 'QGTaggerPFCHS' + label,  
             QGTagger.clone( srcJets   = 'patJetsAK4PFCHSLeg' + label ,jetsLabel = 'QGL_AK4PFchs', vertexIndex = cms.uint32(vertexIndex),
                             srcVertexCollection = 'offlineSlimmedPrimaryVertices'))
  
    from RecoJets.JetProducers.PileupJetIDParams_cfi import full_80x_chs
    from RecoJets.JetProducers.PileupJetIDParams_cfi import full_81x_chs
    if os.environ["CMSSW_VERSION"].count("CMSSW_9"):
        pujidparam = full_81x_chs
    if os.environ["CMSSW_VERSION"].count("CMSSW_8"):
        pujidparam = full_80x_chs
    
    flashggJets = cms.EDProducer('FlashggJetProducer',
                                 DiPhotonTag = cms.InputTag('flashggDiPhotons'),
                                 VertexTag   = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                 JetTag      = cms.InputTag('patJetsAK4PFCHSLeg' + label),
                                 VertexCandidateMapTag = cms.InputTag("flashggVertexMapForCHS"),
                                 qgVariablesInputTag   = cms.InputTag('QGTaggerPFCHS'+label, 'qgLikelihood'),
                                 ComputeSimpleRMS = cms.bool(True),
                                 PileupJetIdParameters = full_81x_chs,
                                 rho     = cms.InputTag("fixedGridRhoFastjetAll"),
                                 JetCollectionIndex = cms.uint32(vertexIndex),
                                 Debug = cms.untracked.bool(debug),
                                 DoPuJetID = cms.bool(True),
                                 ComputeRegVars = cms.bool(False),
                                 MinPtForEneSum = cms.double(0.),
                                 MaxEtaForEneSum = cms.double(2.5),
                                 NJetsForEneSum = cms.uint32(0),
                                 MiniAodJetTag = cms.InputTag("slimmedJets")
                                 )
    setattr( process, 'flashggPFCHSJets'+ label, flashggJets)
  
    # randomize Jets
    from flashgg.MicroAOD.flashggRandomizedJetProducer_cfi import flashggRandomizedJets
    flashggRandomizedPFCHSJets = flashggRandomizedJets.clone()
    flashggRandomizedPFCHSJets.src = "flashggPFCHSJets" + label
    setattr(process.RandomNumberGeneratorService, 'flashggRandomizedPFCHSJets' + label, cms.PSet(initialSeed = cms.untracked.uint32(36423784 + int(label))))
    setattr(process, 'flashggRandomizedPFCHSJets' + label, flashggRandomizedPFCHSJets )
  
    flashggSelectedJets = cms.EDFilter("FLASHggJetSelector",
                                       src = cms.InputTag( 'flashggRandomizedPFCHSJets' + label ),
                                       cut = cms.string("pt > 15.")
    )
    setattr( process, 'flashggSelectedPFCHSJets'+label, flashggSelectedJets )

def buildFinalJetsTask (process, isData):

    for vtx in range(0,maxJetCollections):
        addFlashggPFCHSJets (process = process,
                             isData=isData,
                             vertexIndex =vtx,
                            #doQGTagging = True,
                             label = '' + str(vtx),
                             )

    if isData:
        for vtx in range(0,maxJetCollections):
            delattr(process,"patJetGenJetMatchAK4PFCHSLeg%i"%vtx)
            delattr(process,"patJetFlavourAssociationAK4PFCHSLeg%i"%vtx)
            delattr(process,"patJetPartons%i"%vtx)
            delattr(process,"patJetPartonMatchAK4PFCHSLeg%i"%vtx)
    
    JetCollectionVInputTag      = cms.VInputTag()
    for i in range(0,maxJetCollections):
        JetCollectionVInputTag.append(cms.InputTag('flashggSelectedPFCHSJets' + str(i)))
    
    flashggFinalJets = cms.EDProducer("FlashggVectorVectorJetCollector",
                                      inputTagJets = JetCollectionVInputTag
    )
    
    setattr(process, 'flashggFinalJets', flashggFinalJets)

    #setattr( process, 'jetTask', cms.Task() )
    #getattr( process, 'jetTask', cms.Task() ).add(*[getattr(process,prod) for prod in process.producers_()])
    #getattr( process, 'jetTask', cms.Task() ).add(*[getattr(process,filt) for filt in process.filters_()])



