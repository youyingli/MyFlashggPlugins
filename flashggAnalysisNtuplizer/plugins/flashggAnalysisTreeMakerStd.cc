
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/DataFormats/interface/Met.h"
#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "flashgg/Taggers/interface/LeptonSelection.h"

#include "MyFlashggPlugins/flashggAnalysisNtuplizer/interface/DataFormats.h"
#include "TLorentzVector.h"
#include "TTree.h"

// **********************************************************************

using namespace std;
using namespace edm;
using namespace flashgg;


// **********************************************************************

class flashggAnalysisTreeMakerStd : public edm::EDAnalyzer
{

    public:
        explicit flashggAnalysisTreeMakerStd( const edm::ParameterSet & );
        ~flashggAnalysisTreeMakerStd();
    
        static void fillDescriptions( edm::ConfigurationDescriptions &descriptions );
    
    
    private:
    
        typedef std::vector<edm::Handle<edm::View<flashgg::Jet> > > JetCollectionVector;
        edm::Service<TFileService> fs_;
        TTree* tree;
    
        virtual void beginJob() override;
        virtual void analyze( const edm::Event &, const edm::EventSetup & ) override;
        virtual void endJob() override;
    
        edm::EDGetTokenT<edm::View<flashgg::DiPhotonCandidate>> diphotonToken_;
        std::vector<edm::InputTag> inputTagJets_;
        std::vector<edm::EDGetTokenT<View<flashgg::Jet>>>       tokenJets_;
        edm::EDGetTokenT<edm::View<flashgg::Electron>>          electronToken_;
        edm::EDGetTokenT<edm::View<flashgg::Muon>>              muonToken_;
        edm::EDGetTokenT<edm::View<flashgg::Met>>               metToken_;
    
        edm::EDGetTokenT<edm::View<reco::Vertex>>               vertexToken_;
        edm::EDGetTokenT<reco::BeamSpot>                        beamSpotToken_;
        edm::EDGetTokenT<double>                                rhoTaken_;
        edm::EDGetTokenT<View<reco::GenParticle>>               genParticleToken_;
        edm::EDGetTokenT<GenEventInfoProduct>                   genEventInfoToken_;
        edm::EDGetTokenT<View<PileupSummaryInfo>>               pileUpToken_;
        edm::EDGetTokenT<edm::TriggerResults>                   triggerToken_;
        edm::EDGetTokenT<edm::TriggerResults>                   mettriggerToken_;
        std::string pathName_;
        bool isMiniAOD_;
        bool doHTXS_;
    
        edm::EDGetTokenT<int> stage0catToken_, stage1catToken_, njetsToken_;
        edm::EDGetTokenT<HTXS::HiggsClassification> newHTXSToken_;
        edm::EDGetTokenT<float> pTHToken_,pTVToken_;
    
        flashggAnalysisTreeFormatStd dataformat;

};

// ******************************************************************************************


//
// constructors and destructor
//
flashggAnalysisTreeMakerStd::flashggAnalysisTreeMakerStd( const edm::ParameterSet &iConfig ):
    diphotonToken_        ( consumes< View<flashgg::DiPhotonCandidate> >     ( iConfig.getParameter<InputTag> ( "DiPhotonTag"    ) ) ),
    inputTagJets_         ( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" )                                      ),
    electronToken_        ( consumes< View<flashgg::Electron> >              ( iConfig.getParameter<InputTag> ( "ElectronTag"    ) ) ),
    muonToken_            ( consumes< View<flashgg::Muon> >                  ( iConfig.getParameter<InputTag> ( "MuonTag"        ) ) ),
    metToken_             ( consumes< View<flashgg::Met> >                   ( iConfig.getParameter<InputTag> ( "MetTag"         ) ) ),
    vertexToken_          ( consumes< View<reco::Vertex> >                   ( iConfig.getParameter<InputTag> ( "VertexTag"      ) ) ),
    beamSpotToken_        ( consumes< reco::BeamSpot >                       ( iConfig.getParameter<InputTag> ( "BeamSpotTag"    ) ) ),
    rhoTaken_             ( consumes< double >                               ( iConfig.getParameter<InputTag> ( "RhoTag"         ) ) ),
    genParticleToken_     ( consumes< View<reco::GenParticle> >              ( iConfig.getParameter<InputTag> ( "GenParticleTag" ) ) ),
    genEventInfoToken_    ( consumes< GenEventInfoProduct >                  ( iConfig.getParameter<InputTag> ( "GenEventInfo"   ) ) ),
    pileUpToken_          ( consumes< View<PileupSummaryInfo > >             ( iConfig.getParameter<InputTag> ( "PileUpTag"      ) ) ),
    triggerToken_         ( consumes< edm::TriggerResults >                  ( iConfig.getParameter<InputTag> ( "TriggerTag"     ) ) ),
    mettriggerToken_      ( consumes< edm::TriggerResults >                  ( iConfig.getParameter<InputTag> ( "MetTriggerTag"  ) ) )
{
    pathName_   = iConfig.getParameter<std::string>( "pathName" ) ;
    isMiniAOD_  = iConfig.getParameter<bool>( "isMiniAOD" ) ;
    doHTXS_     = iConfig.getParameter<bool>( "doHTXS" ) ;

    for (unsigned i = 0 ; i < inputTagJets_.size() ; i++) {
        auto token = consumes<View<flashgg::Jet> >(inputTagJets_[i]);
        tokenJets_.push_back(token);
    }

    ParameterSet HTXSps = iConfig.getParameterSet( "HTXSTags" );
    stage0catToken_ = consumes<int>( HTXSps.getParameter<InputTag>("stage0cat") );
    stage1catToken_ = consumes<int>( HTXSps.getParameter<InputTag>("stage1cat") );
    njetsToken_     = consumes<int>( HTXSps.getParameter<InputTag>("njets") );
    pTHToken_       = consumes<float>( HTXSps.getParameter<InputTag>("pTH") );
    pTVToken_       = consumes<float>( HTXSps.getParameter<InputTag>("pTV") );
    newHTXSToken_   = consumes<HTXS::HiggsClassification>( HTXSps.getParameter<InputTag>("ClassificationObj") );

}

flashggAnalysisTreeMakerStd::~flashggAnalysisTreeMakerStd()
{
}

void
flashggAnalysisTreeMakerStd::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup )
{

    // Access edm objects
    // ---------------------------------------------------------------------------------------------------------
    Handle< View<flashgg::DiPhotonCandidate> >   diphotons;
    JetCollectionVector Jets( inputTagJets_.size() );
    Handle< View<flashgg::Electron> >            electrons;
    Handle< View<flashgg::Muon> >                muons;
    Handle< View<flashgg::Met> >                 met;
    Handle< View<reco::Vertex> >                 primaryVertices;
    Handle< reco::BeamSpot >                     recoBeamSpotHandle;
    Handle< double >                             rho;
    Handle< View<reco::GenParticle> >            genParticles;
    Handle< GenEventInfoProduct >                genEventInfo;
    Handle< View< PileupSummaryInfo> >           pileupInfo;
    Handle< edm::TriggerResults >                triggerHandle;
    Handle< edm::TriggerResults >                mettriggerHandle;

    iEvent.getByToken( diphotonToken_     ,     diphotons           );
    for( size_t j = 0; j < inputTagJets_.size(); ++j ) iEvent.getByToken( tokenJets_[j], Jets[j] );
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

    Handle<int> stage0cat, stage1cat, njets;
    Handle<float> pTH, pTV;
    iEvent.getByToken(stage0catToken_, stage0cat);
    iEvent.getByToken(stage1catToken_, stage1cat);
    iEvent.getByToken(njetsToken_, njets);
    iEvent.getByToken(pTHToken_, pTH);
    iEvent.getByToken(pTVToken_, pTV);

    Handle<HTXS::HiggsClassification> htxsClassification;
    iEvent.getByToken(newHTXSToken_, htxsClassification);

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
    auto passMETFilter = [ &mettriggerHandle, &mettriggername ] ( const std::string& triggername ) {
         const unsigned index = mettriggername.triggerIndex( triggername );
         return mettriggerHandle->accept( index ) && mettriggerHandle->wasrun( index ) && !mettriggerHandle->error( index );
    };
                                                                                                                             
    dataformat.Flag_HBHENoiseFilter                    = passMETFilter("Flag_HBHENoiseFilter");
    dataformat.Flag_HBHENoiseIsoFilter                 = passMETFilter("Flag_HBHENoiseIsoFilter");
    dataformat.Flag_EcalDeadCellTriggerPrimitiveFilter = passMETFilter("Flag_EcalDeadCellTriggerPrimitiveFilter");
    dataformat.Flag_goodVertices                       = passMETFilter("Flag_goodVertices");
    dataformat.Flag_globalSuperTightHalo2016Filter     = passMETFilter("Flag_globalSuperTightHalo2016Filter");
    dataformat.Flag_BadPFMuonFilter                    = passMETFilter("Flag_BadPFMuonFilter");
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
            dataformat.elecs_fggPhoVeto          .emplace_back( phoVeto( it_elec, diphoPtr, 0., 0.4, 0. ) );

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

    //Athough choosing leading diphoton information and storing them and associated ones, 
    //but still keep global information and gen information for completeness.
    tree->Fill();

}


void
flashggAnalysisTreeMakerStd::beginJob()
{
    tree = fs_->make<TTree>( "flashggStdTree", "flashgg Standard Analysis Tree" );
    dataformat.RegisterTree( tree );
}

void
flashggAnalysisTreeMakerStd::endJob()
{
}

void
flashggAnalysisTreeMakerStd::fillDescriptions( edm::ConfigurationDescriptions &descriptions )
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault( desc );
}

DEFINE_FWK_MODULE( flashggAnalysisTreeMakerStd );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
