#include "DataFormats/Common/interface/Ptr.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "flashgg/Taggers/interface/LeptonSelection.h"

#include "MyFlashggPlugins/flashggAnalysisNtuplizer/interface/flashggAnaTreeMakerWithSyst.h"

using namespace std;
using namespace edm;
using namespace flashgg;

flashggAnaTreeMakerWithSyst::flashggAnaTreeMakerWithSyst( const edm::InputTag &diphoton, const edm::ParameterSet &iConfig, edm::ConsumesCollector&& iC):
    diphotonToken_        ( iC.consumes< View<flashgg::DiPhotonCandidate> >     ( diphoton                                            ) ),
    inputTagJets_         ( iConfig.getParameter<vector<edm::InputTag> >( "inputTagJets" )                                              ),
    electronToken_        ( iC.consumes< View<flashgg::Electron> >              ( iConfig.getParameter<InputTag> ( "ElectronTag"    ) ) ),
    muonToken_            ( iC.consumes< View<flashgg::Muon> >                  ( iConfig.getParameter<InputTag> ( "MuonTag"        ) ) ),
    metToken_             ( iC.consumes< View<flashgg::Met> >                   ( iConfig.getParameter<InputTag> ( "MetTag"         ) ) ),
    vertexToken_          ( iC.consumes< View<reco::Vertex> >                   ( iConfig.getParameter<InputTag> ( "VertexTag"      ) ) ),
    beamSpotToken_        ( iC.consumes< reco::BeamSpot >                       ( iConfig.getParameter<InputTag> ( "BeamSpotTag"    ) ) ),
    rhoTaken_             ( iC.consumes< double >                               ( iConfig.getParameter<InputTag> ( "RhoTag"         ) ) ),
    genParticleToken_     ( iC.consumes< View<reco::GenParticle> >              ( iConfig.getParameter<InputTag> ( "GenParticleTag" ) ) ),
    genEventInfoToken_    ( iC.consumes< GenEventInfoProduct >                  ( iConfig.getParameter<InputTag> ( "GenEventInfo"   ) ) ),
    pileUpToken_          ( iC.consumes< View<PileupSummaryInfo > >             ( iConfig.getParameter<InputTag> ( "PileUpTag"      ) ) ),
    triggerToken_         ( iC.consumes< edm::TriggerResults >                  ( iConfig.getParameter<InputTag> ( "TriggerTag"     ) ) ),
    mettriggerToken_      ( iC.consumes< edm::TriggerResults >                  ( iConfig.getParameter<InputTag> ( "MetTriggerTag"  ) ) )
{
    pathName_   = iConfig.getParameter<std::string>( "pathName" ) ;
    isMiniAOD_  = iConfig.getParameter<bool>( "isMiniAOD" ) ;
    doHTXS_     = iConfig.getParameter<bool>( "doHTXS" ) ;

    for (unsigned i = 0 ; i < inputTagJets_.size() ; i++) {
        auto token = iC.consumes<View<flashgg::Jet> >(inputTagJets_[i]);
        tokenJets_.push_back(token);
    }

    ParameterSet HTXSps = iConfig.getParameterSet( "HTXSTags" );
    newHTXSToken_ = iC.consumes<HTXS::HiggsClassification>( HTXSps.getParameter<InputTag>("ClassificationObj") );

}

flashggAnaTreeMakerWithSyst::~flashggAnaTreeMakerWithSyst()
{
}

void
flashggAnaTreeMakerWithSyst::RegisterTree( TTree* tree )
{
    dataformat.RegisterTree( tree );
}

void
flashggAnaTreeMakerWithSyst::Analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup )
{

    // Access edm objects
    // ---------------------------------------------------------------------------------------------------------
    JetCollectionVector Jets( inputTagJets_.size() );

    iEvent.getByToken( diphotonToken_     ,     diphotons           );
    for( size_t j = 0; j < inputTagJets_.size(); ++j ) 
        iEvent.getByToken( tokenJets_[j], Jets[j] );
    iEvent.getByToken( electronToken_     ,     electrons           );
    iEvent.getByToken( muonToken_         ,     muons               );
    iEvent.getByToken( metToken_          ,     met                 );
    iEvent.getByToken( vertexToken_       ,     primaryVertices     );
    iEvent.getByToken( beamSpotToken_     ,     recoBeamSpotHandle  );
    iEvent.getByToken( rhoTaken_          ,     rho                 );
    iEvent.getByToken( genParticleToken_  ,     genParticles        );
    iEvent.getByToken( genEventInfoToken_ ,     genEventInfo        );
    iEvent.getByToken( pileUpToken_       ,     pileupInfo          );
    iEvent.getByToken( triggerToken_      ,     triggerHandle       );
    iEvent.getByToken( mettriggerToken_   ,     mettriggerHandle    );
    iEvent.getByToken( newHTXSToken_      ,     htxsClassification  );

    // dataformat Initialzation
    // ---------------------------------------------------------------------------------------------------------
    dataformat.Initialzation();

    // Global information
    // ---------------------------------------------------------------------------------------------------------
    dataformat.Rho     = *rho;
    dataformat.PVz     = primaryVertices->ptrAt(0)->z();
    dataformat.NVtx    = primaryVertices->size();
    if( recoBeamSpotHandle.isValid() ) {
        dataformat.BSsigmaz = recoBeamSpotHandle->sigmaZ();
    }

    if (iEvent.isRealData()) {
        dataformat.passTrigger = true;
    } else {
        const edm::TriggerNames &triggerNames = iEvent.triggerNames( *triggerHandle );
                                                                                                             
        map<string, int> triggerIndices;
        for( unsigned int i = 0; i < triggerNames.triggerNames().size(); i++ ) {
            std::string trimmedName = HLTConfigProvider::removeVersion( triggerNames.triggerName( i ) );
            triggerIndices.emplace(trimmedName, triggerNames.triggerIndex( triggerNames.triggerName( i ) ));
        }
                                                                                                             
        for (const auto& it : triggerIndices) {
            if (triggerHandle->accept(it.second)) {
                if (it.first == pathName_) {dataformat.passTrigger = true; break;}
            }
        }
    }

    const edm::TriggerNames& mettriggername = iEvent.triggerNames( *mettriggerHandle );
    auto passMETFilter = [ this, &mettriggername ] ( const std::string& triggername ) {
         const unsigned index = mettriggername.triggerIndex( triggername );
         return    this->mettriggerHandle->accept( index ) 
                && this->mettriggerHandle->wasrun( index ) 
                && !this->mettriggerHandle->error( index );
    };
                                                                                                                             
    dataformat.Flag_HBHENoiseFilter                    = passMETFilter("Flag_HBHENoiseFilter");
    dataformat.Flag_HBHENoiseIsoFilter                 = passMETFilter("Flag_HBHENoiseIsoFilter");
    dataformat.Flag_EcalDeadCellTriggerPrimitiveFilter = passMETFilter("Flag_EcalDeadCellTriggerPrimitiveFilter");
    dataformat.Flag_goodVertices                       = passMETFilter("Flag_goodVertices");
    dataformat.Flag_globalSuperTightHalo2016Filter     = passMETFilter("Flag_globalSuperTightHalo2016Filter");
    dataformat.Flag_BadPFMuonFilter                    = passMETFilter("Flag_BadPFMuonFilter");
    dataformat.Flag_BadChargedCandidateFilter          = passMETFilter("Flag_BadChargedCandidateFilter");
    dataformat.Flag_ecalBadCalibFilter                 = passMETFilter("Flag_ecalBadCalibFilter");
    dataformat.Flag_eeBadScFilter                      = iEvent.isRealData() ? passMETFilter("Flag_eeBadScFilter") : true;

    // Choose leading diphoton information and store them and associated ones
    const std::vector<edm::Ptr<flashgg::DiPhotonCandidate> > diphotonPtrs = diphotons->ptrs();
    if (diphotonPtrs.size() > 0) {

        // DiPhoton information 
        // ---------------------------------------------------------------------------------------------------------
        const Ptr<flashgg::DiPhotonCandidate> diphoPtr = diphotonPtrs[0];
        dataformat.dipho_mass                 = diphoPtr->mass();
        dataformat.dipho_pt                   = diphoPtr->pt();
        dataformat.dipho_leadPt               = diphoPtr->leadingPhoton()->pt();
        dataformat.dipho_leadEta              = diphoPtr->leadingPhoton()->eta();
        dataformat.dipho_leadPhi              = diphoPtr->leadingPhoton()->phi();
        dataformat.dipho_leadE                = diphoPtr->leadingPhoton()->energy();
        dataformat.dipho_leadsigEOverE        = diphoPtr->leadingPhoton()->sigEOverE();
        dataformat.dipho_leadR9               = diphoPtr->leadingPhoton()->full5x5_r9();
        dataformat.dipho_leadsieie            = diphoPtr->leadingPhoton()->full5x5_sigmaIetaIeta();
        dataformat.dipho_leadhoe              = diphoPtr->leadingPhoton()->hadronicOverEm();
        dataformat.dipho_leadIDMVA            = diphoPtr->leadingView()->phoIdMvaWrtChosenVtx();
        dataformat.dipho_leadhasPixelSeed     = diphoPtr->leadingPhoton()->hasPixelSeed();
        dataformat.dipho_leadGenMatch         = diphoPtr->leadingPhoton()->hasMatchedGenPhoton();
        dataformat.dipho_leadGenMatchType     = diphoPtr->leadingPhoton()->genMatchType();//enum mcMatch_t { kUnkown = 0, kPrompt, kFake  };
        dataformat.dipho_subleadPt            = diphoPtr->subLeadingPhoton()->pt();
        dataformat.dipho_subleadEta           = diphoPtr->subLeadingPhoton()->eta();
        dataformat.dipho_subleadPhi           = diphoPtr->subLeadingPhoton()->phi();
        dataformat.dipho_subleadE             = diphoPtr->subLeadingPhoton()->energy();
        dataformat.dipho_subleadsigEOverE     = diphoPtr->subLeadingPhoton()->sigEOverE();
        dataformat.dipho_subleadR9            = diphoPtr->subLeadingPhoton()->full5x5_r9();
        dataformat.dipho_subleadsieie         = diphoPtr->subLeadingPhoton()->full5x5_sigmaIetaIeta();
        dataformat.dipho_subleadhoe           = diphoPtr->subLeadingPhoton()->hadronicOverEm();
        dataformat.dipho_subleadIDMVA         = diphoPtr->subLeadingView()->phoIdMvaWrtChosenVtx();
        dataformat.dipho_subleadhasPixelSeed  = diphoPtr->subLeadingPhoton()->hasPixelSeed();
        dataformat.dipho_subleadGenMatch      = diphoPtr->subLeadingPhoton()->hasMatchedGenPhoton();
        dataformat.dipho_subleadGenMatchType  = diphoPtr->subLeadingPhoton()->genMatchType();
        dataformat.dipho_SelectedVz           = diphoPtr->vtx()->position().z();
        dataformat.dipho_GenVz                = diphoPtr->genPV().z();

        // Electron information
        // ---------------------------------------------------------------------------------------------------------
        int Nelecs = 0;
        for ( const auto& it_elec : electrons->ptrs() ) {

            double elecEta = fabs( it_elec->superCluster()->eta() );
            if ( elecEta > 2.5 || ( elecEta > 1.4442 && elecEta < 1.566 ) ) continue;

            dataformat.elecs_Charge              .emplace_back( it_elec->charge() );
            dataformat.elecs_Pt                  .emplace_back( it_elec->pt() );
            dataformat.elecs_Eta                 .emplace_back( it_elec->eta() );
            dataformat.elecs_Phi                 .emplace_back( it_elec->phi() );
            dataformat.elecs_Energy              .emplace_back( it_elec->energy() );
            dataformat.elecs_EtaSC               .emplace_back( it_elec->superCluster()->eta() );
            dataformat.elecs_PhiSC               .emplace_back( it_elec->superCluster()->phi() );
            dataformat.elecs_EGMCutBasedIDVeto   .emplace_back( it_elec->passVetoId() );
            dataformat.elecs_EGMCutBasedIDLoose  .emplace_back( it_elec->passLooseId() );
            dataformat.elecs_EGMCutBasedIDMedium .emplace_back( it_elec->passMediumId() );
            dataformat.elecs_EGMCutBasedIDTight  .emplace_back( it_elec->passTightId() );
            dataformat.elecs_fggPhoVeto          .emplace_back( phoVeto( it_elec, diphoPtr, 0.2, 0.35, 0.0 ) );
            dataformat.elecs_tmpPhoVeto          .emplace_back( phoVeto( it_elec, diphoPtr, 0.4, 0.4, 0.0 ) );

            if ( isMiniAOD_ ) {
                dataformat.elecs_GsfTrackDz     .emplace_back( it_elec->gsfTrack()->dz( diphoPtr->vtx()->position() ) );
                dataformat.elecs_GsfTrackDxy    .emplace_back( it_elec->gsfTrack()->dxy( diphoPtr->vtx()->position() ) );
                if ( !iEvent.isRealData() ) {
                    const reco::GenParticle* gen = it_elec->genLepton();
                    if ( gen != nullptr ) {
                        dataformat.elecs_GenMatch .emplace_back( true );
                        dataformat.elecs_GenPdgID .emplace_back( gen->pdgId() );
                        dataformat.elecs_GenPt    .emplace_back( gen->pt() );
                        dataformat.elecs_GenEta   .emplace_back( gen->eta() );
                        dataformat.elecs_GenPhi   .emplace_back( gen->phi() );
                    } else {
                        dataformat.elecs_GenMatch .emplace_back( false );
                        dataformat.elecs_GenPdgID .emplace_back( 0 );
                        dataformat.elecs_GenPt    .emplace_back( -999. );
                        dataformat.elecs_GenEta   .emplace_back( -999. );
                        dataformat.elecs_GenPhi   .emplace_back( -999. );
                    }
                }
            }

            Nelecs++;
        }
        dataformat.elecs_size = Nelecs;

        // Muon information
        // ---------------------------------------------------------------------------------------------------------
        double Nmuons = 0;
        for ( const auto& it_muon : muons->ptrs() ) {
            if ( !it_muon->passed(reco::Muon::CutBasedIdLoose) ) continue;
            if ( fabs( it_muon->eta() ) > 2.4 ) continue;

            float dRPho1Muon = reco::deltaR( it_muon->eta(), it_muon->phi(), dataformat.dipho_leadEta, dataformat.dipho_leadPhi );
            float dRPho2Muon = reco::deltaR( it_muon->eta(), it_muon->phi(), dataformat.dipho_subleadEta, dataformat.dipho_subleadPhi );
            if ( dRPho1Muon < 0.2 || dRPho2Muon < 0.2 ) continue;

            dataformat.muons_Charge                  .emplace_back( it_muon->charge() );
            dataformat.muons_MuonType                .emplace_back( it_muon->type() );
            dataformat.muons_Pt                      .emplace_back( it_muon->pt() );
            dataformat.muons_Eta                     .emplace_back( it_muon->eta() );
            dataformat.muons_Phi                     .emplace_back( it_muon->phi() );
            dataformat.muons_Energy                  .emplace_back( it_muon->energy() );
            dataformat.muons_BestTrackDz             .emplace_back( it_muon->muonBestTrack()->dz( diphoPtr->vtx()->position() ) );
            dataformat.muons_BestTrackDxy            .emplace_back( it_muon->muonBestTrack()->dxy( diphoPtr->vtx()->position() ) );
            dataformat.muons_CutBasedIdMedium        .emplace_back( it_muon->passed(reco::Muon::CutBasedIdMedium) );
            dataformat.muons_CutBasedIdTight         .emplace_back( muon::isTightMuon( *it_muon, *(diphoPtr->vtx()) ) );
            dataformat.muons_PFIsoDeltaBetaCorrR04   .emplace_back( it_muon->fggPFIsoSumRelR04() );
            dataformat.muons_TrackerBasedIsoR03      .emplace_back( it_muon->fggTrkIsoSumRelR03() );

            if ( isMiniAOD_ ) {
                if ( !iEvent.isRealData() ) {
                    const reco::GenParticle* gen = it_muon->genLepton();
                    if ( gen != nullptr ) {
                        dataformat.muons_GenMatch .emplace_back( true );
                        dataformat.muons_GenPdgID .emplace_back( gen->pdgId() );
                        dataformat.muons_GenPt    .emplace_back( gen->pt() );
                        dataformat.muons_GenEta   .emplace_back( gen->eta() );
                        dataformat.muons_GenPhi   .emplace_back( gen->phi() );
                    } else {
                        dataformat.muons_GenMatch .emplace_back( false );
                        dataformat.muons_GenPdgID .emplace_back( 0 );
                        dataformat.muons_GenPt    .emplace_back( -999. );
                        dataformat.muons_GenEta   .emplace_back( -999. );
                        dataformat.muons_GenPhi   .emplace_back( -999. );
                    }
                }
            }

            Nmuons++;
        }
        dataformat.muons_size = Nmuons;

        // Jet information
        // ---------------------------------------------------------------------------------------------------------
        int Njets = 0;
        unsigned int jetCollectionIndex = diphoPtr->jetCollectionIndex();
        for ( const auto& it_jet : Jets[jetCollectionIndex]->ptrs() ) {
            if ( !it_jet->passesJetID(flashgg::Tight2017) ) continue;
            if ( fabs( it_jet->eta() ) > 4.7 ) { continue; }

            float dRPho1Jet = reco::deltaR( it_jet->eta(), it_jet->phi(), dataformat.dipho_leadEta, dataformat.dipho_leadPhi );
            float dRPho2Jet = reco::deltaR( it_jet->eta(), it_jet->phi(), dataformat.dipho_subleadEta, dataformat.dipho_subleadPhi );
            if ( dRPho1Jet < 0.4 || dRPho2Jet < 0.4 ) continue;

            dataformat.jets_Pt                            .emplace_back( it_jet->pt() );
            dataformat.jets_Eta                           .emplace_back( it_jet->eta() );
            dataformat.jets_Phi                           .emplace_back( it_jet->phi() );
            dataformat.jets_Mass                          .emplace_back( it_jet->mass() );
            dataformat.jets_Energy                        .emplace_back( it_jet->energy() );
            dataformat.jets_PtRaw                         .emplace_back( it_jet->correctedJet( "Uncorrected" ).pt() );
            dataformat.jets_QGL                           .emplace_back( it_jet->QGL() );
            dataformat.jets_RMS                           .emplace_back( it_jet->rms() );
            dataformat.jets_puJetIdMVA                    .emplace_back( it_jet->puJetIdMVA() );
            dataformat.jets_GenJetMatch                   .emplace_back( it_jet->hasGenMatch() );
            dataformat.jets_pfCombinedInclusiveSecondaryVertexV2BJetTags        
                                                          .emplace_back( it_jet->bDiscriminator( "pfCombinedInclusiveSecondaryVertexV2BJetTags" ) );
            dataformat.jets_pfCombinedMVAV2BJetTags       .emplace_back( it_jet->bDiscriminator( "pfCombinedMVAV2BJetTags" ) );
            dataformat.jets_pfDeepCSVJetTags_probb        .emplace_back( it_jet->bDiscriminator( "pfDeepCSVJetTags:probb" ) );
            dataformat.jets_pfDeepCSVJetTags_probbb       .emplace_back( it_jet->bDiscriminator( "pfDeepCSVJetTags:probbb" ) );
            dataformat.jets_pfDeepCSVJetTags_probc        .emplace_back( it_jet->bDiscriminator( "pfDeepCSVJetTags:probc" ) );
            dataformat.jets_pfDeepCSVJetTags_probudsg     .emplace_back( it_jet->bDiscriminator( "pfDeepCSVJetTags:probudsg" ) );

            if ( isMiniAOD_ && !iEvent.isRealData() ) {
                const reco::GenParticle* parton = it_jet->genParton();
                if ( parton != nullptr ) {
                    dataformat.jets_GenPartonMatch    .emplace_back( true );
                    dataformat.jets_GenPt             .emplace_back( parton->pt() );
                    dataformat.jets_GenEta            .emplace_back( parton->eta() );
                    dataformat.jets_GenPhi            .emplace_back( parton->phi() );
                    dataformat.jets_GenPdgID          .emplace_back( parton->pdgId() );
                    dataformat.jets_GenFlavor         .emplace_back( it_jet->partonFlavour() );
                    dataformat.jets_GenHadronFlavor   .emplace_back( it_jet->hadronFlavour() );
                } else {
                    dataformat.jets_GenPartonMatch    .emplace_back( false );
                    dataformat.jets_GenPt             .emplace_back( -999. );
                    dataformat.jets_GenEta            .emplace_back( -999. );
                    dataformat.jets_GenPhi            .emplace_back( -999. );
                    dataformat.jets_GenPdgID          .emplace_back( 0 );
                    dataformat.jets_GenFlavor         .emplace_back( 0 );
                    dataformat.jets_GenHadronFlavor   .emplace_back( 0 );
                }
            }

            Njets++;
        } //jet loop
        dataformat.jets_size = Njets;

        // Met information
        // ---------------------------------------------------------------------------------------------------------
        edm::Ptr<flashgg::Met> theMet = met->ptrAt( 0 );
        dataformat.met_Pt     = theMet->pt();
        dataformat.met_Phi    = theMet->phi();
        dataformat.met_Px     = theMet->px();
        dataformat.met_Py     = theMet->py();
        dataformat.met_SumET  = theMet->sumEt();

    }

    // Gen information
    // ---------------------------------------------------------------------------------------------------------
    if (!iEvent.isRealData()) {
        for( unsigned int PVI = 0; PVI < pileupInfo->size(); ++PVI ) {
            Int_t pu_bunchcrossing = pileupInfo->ptrAt( PVI )->getBunchCrossing();
            if( pu_bunchcrossing == 0 ) { dataformat.NPu = pileupInfo->ptrAt( PVI )->getTrueNumInteractions(); break; }
        }

        dataformat.genweight     = genEventInfo->weight();

        int NGenParticles = 0;
        const std::vector<edm::Ptr<reco::GenParticle> > genParticlesPtrs = genParticles->ptrs();
        for (const auto& it_gen : genParticles->ptrs()) {
            if ( abs(it_gen->pdgId()) > 25) continue;
            if ( it_gen->status() > 30 ) continue;
            if ( it_gen->pdgId() == 22      && it_gen->status() == 1 && !(it_gen->pt() > 10 || it_gen->isPromptFinalState())) continue;
            if ( abs(it_gen->pdgId()) == 11 && it_gen->status() == 1 && !(it_gen->pt() > 3  || it_gen->isPromptFinalState())) continue;

            dataformat.GenParticles_Pt     .emplace_back( it_gen->pt() );
            dataformat.GenParticles_Eta    .emplace_back( it_gen->eta() );
            dataformat.GenParticles_Phi    .emplace_back( it_gen->phi() );
            dataformat.GenParticles_Mass   .emplace_back( it_gen->mass() );
            dataformat.GenParticles_PdgID  .emplace_back( it_gen->pdgId() );
            dataformat.GenParticles_Status .emplace_back( it_gen->status() );
            dataformat.GenParticles_nMo    .emplace_back( it_gen->numberOfMothers() );
            dataformat.GenParticles_nDa    .emplace_back( it_gen->numberOfDaughters() );

            NGenParticles++;
        }

        dataformat.GenParticles_size = NGenParticles;

        if(doHTXS_) {
            dataformat.HTXSstage0cat = htxsClassification->stage0_cat;
            dataformat.HTXSstage1cat = htxsClassification->stage1_cat_pTjet30GeV;
            dataformat.HTXSnjets     = htxsClassification->jets30.size();
            dataformat.HTXSpTH       = htxsClassification->p4decay_higgs.pt();
            dataformat.HTXSpTV       = htxsClassification->p4decay_V.pt();
        }

    }

}

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
