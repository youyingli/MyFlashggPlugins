import FWCore.ParameterSet.Config as cms
from flashgg.MicroAOD.flashggJets_cfi import maxJetCollections

flashggPhotonJet = cms.EDProducer('FlashggPhotonJetProducer',
                                  JetTag=cms.InputTag('slimmedJets'),
                                  PhotonTag=cms.InputTag('flashggPhotons'),
                                  VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                  VertexSelectorName=cms.string("FlashggLegacyVertexSelector"),
                                  VertexCandidateMapTag=cms.InputTag("flashggVertexMapUnique"),
                                  ConversionTag=cms.InputTag("reducedEgamma","reducedConversions"),
                                  ConversionTagSingleLeg=cms.InputTag("reducedEgamma","reducedSingleLegConversions"),
                                  rhoTag = cms.InputTag( "fixedGridRhoFastjetAll" ),
                                  beamSpotTag = cms.InputTag( "offlineBeamSpot" ),

                                  ##Parameters for Legacy Vertex Selector                                                
                                  #vtxId 2012
                                  #vertexIdMVAweightfile = cms.FileInPath("flashgg/MicroAOD/data/2013FinalPaper_VertexID_BDTCat_conversions.weights.xml"),
                                  #vxtProb2012
                                  #vertexProbMVAweightfile = cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_BDTvtxprob2012.weights.xml"),

                                  #vtxId and vtxProb 2015 no single leg
                                  #useSingleLeg=cms.untracked.bool(False),
                                  #vertexIdMVAweightfile = cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_BDTVtxId_noSL_2015.xml"),
                                  #vertexProbMVAweightfile = cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_BDTVtxProb_noSL_2015.xml"),

                                  #vtxId and vtxProb 2015 with single leg
                                  #vertexIdMVAweightfile = cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_BDTVtxId_SL_2015.xml"),
                                  #vertexProbMVAweightfile = cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_BDTVtxProb_SL_2015.xml"),

                                  #vtxId and vtxProb 2016 with single leg
                                  vertexIdMVAweightfile   = cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_BDTVtxId_SL_2016.xml"),
                                  vertexProbMVAweightfile = cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_BDTVtxProb_SL_2016.xml"),

                                  useSingleLeg=cms.bool(True),
                                  useZerothVertexFromMicro = cms.bool(False),

                                  nVtxSaveInfo = cms.untracked.uint32(999),
                                  trackHighPurity = cms.bool(False),
                                  pureGeomConvMatching    = cms.bool(True),
                                  dRexclude=cms.double(0.05),
                                  #new reso:
                                  sigma1Pix              = cms.double( 0.0125255 ),
                                  sigma1Tib              = cms.double( 0.716301  ),
                                  sigma1Tob              = cms.double( 3.17615   ),
                                  sigma1PixFwd           = cms.double( 0.0581667 ),
                                  sigma1Tid              = cms.double( 0.38521   ),
                                  sigma1Tec              = cms.double( 1.67937   ),
                                  sigma2Pix              = cms.double( 0.0298574 ),
                                  sigma2Tib              = cms.double( 0.414393  ),
                                  sigma2Tob              = cms.double( 1.06805   ),
                                  sigma2PixFwd           = cms.double( 0.180419  ),
                                  sigma2Tid              = cms.double( 0.494722  ),
                                  sigma2Tec              = cms.double( 1.21941   ),
                                  singlelegsigma1Pix     = cms.double( 0.0178107 ),
                                  singlelegsigma1Tib     = cms.double( 1.3188    ),
                                  singlelegsigma1Tob     = cms.double( 2.23662   ),
                                  singlelegsigma1PixFwd  = cms.double( 0.152157  ),
                                  singlelegsigma1Tid     = cms.double( 0.702755  ),
                                  singlelegsigma1Tec     = cms.double( 2.46599   ),
                                  singlelegsigma2Pix     = cms.double( 0.0935307 ),
                                  singlelegsigma2Tib     = cms.double( 0.756568  ),
                                  singlelegsigma2Tob     = cms.double( 0.62143   ),
                                  singlelegsigma2PixFwd  = cms.double( 0.577081  ),
                                  singlelegsigma2Tid     = cms.double( 0.892751  ),
                                  singlelegsigma2Tec     = cms.double( 1.56638   ),
                                  MaxJetCollections = cms.uint32(maxJetCollections),

                                  #Photon jet cut
                                  minJetPt = cms.double(30.),
                                  maxJetEta = cms.double(2.5),
                                  minPhotonPt = cms.double(55.),

                                  #2017 tmp
                                  #Photon loose Id cut https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_94X_and_later
                                  minPhotEBHoE = cms.double(0.04596),
                                  minPhotEEHoE = cms.double(0.059),
                                  minPhotEBsietaieta = cms.double(0.0106),
                                  minPhotEEsietaieta = cms.double(0.0272),
                                  PhotEBChgIsoParams = cms.vdouble(1.694),
                                  PhotEEChgIsoParams = cms.vdouble(2.089),
                                  PhotEBNeuIsoParams = cms.vdouble(24.032, 0.01512, 0.00002259),
                                  PhotEENeuIsoParams = cms.vdouble(19.722, 0.0117,  0.000023),
                                  PhotEBPhoIsoParams = cms.vdouble(2.876, 0.004017),
                                  PhotEEPhoIsoParams = cms.vdouble(4.162, 0.0037),
                                  iphotIsolnAreaValN = cms.int32(7),
                                  photIsolnEAreaVal = cms.vdouble(1.0, 1.479, 2.0, 2.2, 2.3, 2.4, 999999.0),
                                  photIsolnEAreaChgHad = cms.vdouble(0.0112, 0.0108, 0.0106, 0.01002, 0.0098, 0.0089, 0.0087),
                                  photIsolnEAreaNeuHad = cms.vdouble(0.0668, 0.1054, 0.0786, 0.0233, 0.0078, 0.0028, 0.0137),
                                  photIsolnEAreaPhot = cms.vdouble(0.1113, 0.0953, 0.0619, 0.0837, 0.1070, 0.1212, 0.1466)

                                  #2016
                                  #Photon loose Id cut https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_2016_data_for
                                  #minPhotEBHoE = cms.double(0.0597),
                                  #minPhotEEHoE = cms.double(0.0481),
                                  #minPhotEBsietaieta = cms.double(0.01031),
                                  #minPhotEEsietaieta = cms.double(0.03013),
                                  #PhotEBChgIsoParams = cms.vdouble(1.295),
                                  #PhotEEChgIsoParams = cms.vdouble(1.011),
                                  #PhotEBNeuIsoParams = cms.vdouble(10.910, 0.0148, 0.000017),
                                  #PhotEENeuIsoParams = cms.vdouble(5.931, 0.0163, 0.000014),
                                  #PhotEBPhoIsoParams = cms.vdouble(3.630, 0.0047),
                                  #PhotEEPhoIsoParams = cms.vdouble(6.641, 0.0034),
                                  #iphotIsolnAreaValN = cms.int32(7),
                                  #photIsolnEAreaVal = cms.vdouble(1.0, 1.479, 2.0, 2.2, 2.3, 2.4, 999999.0),
                                  #photIsolnEAreaChgHad = cms.vdouble(0.0360, 0.0377, 0.0306, 0.0283, 0.0254, 0.0217, 0.0167),
                                  #photIsolnEAreaNeuHad = cms.vdouble(0.0597, 0.0807, 0.0629, 0.0197, 0.0184, 0.0284, 0.0591),
                                  #photIsolnEAreaPhot = cms.vdouble(0.1210, 0.1107, 0.0699, 0.1056, 0.1457, 0.1719, 0.1998)

                                  )
