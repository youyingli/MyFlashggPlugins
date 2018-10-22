
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
#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "TTree.h"

// **********************************************************************

// define the structures used to create tree branches and fill the trees

struct DiPhotonJetsInfo {

    int nPu                          ;
    float genvz                      ;
    float genweight                  ;

    float rho                        ;
    float pvz                        ;
    float selectedvz                 ;
    int nvertex                      ;
    float BSsigmaz                   ;
    int passTrigger                  ;

    float dipho_mass                 ;
    float dipho_pt                   ;
    float dipho_leadPt               ; 
    float dipho_leadEta              ; 
    float dipho_leadPhi              ; 
    float dipho_leadE                ; 
    float dipho_leadsigEOverE        ; 
    float dipho_leadR9               ; 
    float dipho_leadsieie            ; 
    float dipho_leadhoe              ; 
    float dipho_leadIDMVA            ;
    int   dipho_leadGenMatch         ;
    int   dipho_leadGenMatchType     ;
    float dipho_subleadPt            ; 
    float dipho_subleadEta           ; 
    float dipho_subleadPhi           ; 
    float dipho_subleadE             ; 
    float dipho_subleadsigEOverE     ; 
    float dipho_subleadR9            ; 
    float dipho_subleadsieie         ; 
    float dipho_subleadhoe           ; 
    float dipho_subleadIDMVA         ; 
    int   dipho_subleadGenMatch      ;
    int   dipho_subleadGenMatchType  ;

    int jets_size;
    vector<float> jets_Pt         ;
    vector<float> jets_Eta        ;
    vector<float> jets_Phi        ;
    vector<float> jets_Mass       ;
    vector<float> jets_QGL        ;
    vector<float> jets_RMS        ;
    vector<float> jets_puJetIdMVA ;
    vector<int>   jets_GenMatch   ;

    int HTXSstage0cat ;
    int HTXSstage1cat ;
    int HTXSnjets     ;
    float HTXSpTH     ;
    float HTXSpTV     ; 

};

// **********************************************************************

using namespace std;
using namespace edm;
using namespace flashgg;


// **********************************************************************

class DiPhotonJetsTreeMaker : public edm::EDAnalyzer
{
public:
    explicit DiPhotonJetsTreeMaker( const edm::ParameterSet & );
    ~DiPhotonJetsTreeMaker();

    static void fillDescriptions( edm::ConfigurationDescriptions &descriptions );


private:

    edm::Service<TFileService> fs_;



    virtual void beginJob() override;
    virtual void analyze( const edm::Event &, const edm::EventSetup & ) override;
    virtual void endJob() override;

    void initEventStructure();

    edm::EDGetTokenT<edm::View<flashgg::DiPhotonCandidate> > diphotonToken_;
    std::vector<edm::InputTag> inputTagJets_;
    edm::EDGetTokenT<edm::View<reco::Vertex> >               vertexToken_;

    edm::EDGetTokenT<reco::BeamSpot > beamSpotToken_;
    edm::EDGetTokenT<double> rhoTaken_;
    edm::EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
    edm::EDGetTokenT<GenEventInfoProduct>  genEventInfoToken_;
    edm::EDGetTokenT<View<PileupSummaryInfo>> pileUpToken_;
    edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
    std::string pathName_;
    bool doHTXS_;

    typedef std::vector<edm::Handle<edm::View<flashgg::Jet> > > JetCollectionVector;
    std::vector<edm::EDGetTokenT<View<flashgg::Jet> > > tokenJets_;

    EDGetTokenT<int> stage0catToken_, stage1catToken_, njetsToken_;
    EDGetTokenT<HTXS::HiggsClassification> newHTXSToken_;
    EDGetTokenT<float> pTHToken_,pTVToken_;

    TTree *DiPhotonJetsTree;

    DiPhotonJetsInfo diPhotonJetsInfo;

};

// ******************************************************************************************


//
// constructors and destructor
//
DiPhotonJetsTreeMaker::DiPhotonJetsTreeMaker( const edm::ParameterSet &iConfig ):
    diphotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
    inputTagJets_ ( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) ),
    vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
    beamSpotToken_( consumes<reco::BeamSpot >( iConfig.getParameter<InputTag>( "BeamSpotTag" ) ) ),
    rhoTaken_( consumes<double>( iConfig.getParameter<InputTag>( "RhoTag" ) ) ),
    genParticleToken_( consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag> ( "GenParticleTag" ) ) ),
    genEventInfoToken_( consumes<GenEventInfoProduct>( iConfig.getParameter<InputTag> ( "GenEventInfo" ) ) ),
    pileUpToken_( consumes<View<PileupSummaryInfo> >( iConfig.getParameter<InputTag> ( "PileUpTag" ) ) ),
    triggerToken_( consumes<edm::TriggerResults>( iConfig.getParameter<InputTag> ( "TriggerTag" ) ) )
{
    pathName_ = iConfig.getParameter<std::string>( "pathName" ) ;
    doHTXS_   = iConfig.getParameter<bool>( "doHTXS" ) ;

    for (unsigned i = 0 ; i < inputTagJets_.size() ; i++) {
        auto token = consumes<View<flashgg::Jet> >(inputTagJets_[i]);
        tokenJets_.push_back(token);
    }

    ParameterSet HTXSps = iConfig.getParameterSet( "HTXSTags" );
    stage0catToken_ = consumes<int>( HTXSps.getParameter<InputTag>("stage0cat") );
    stage1catToken_ = consumes<int>( HTXSps.getParameter<InputTag>("stage1cat") );
    njetsToken_ = consumes<int>( HTXSps.getParameter<InputTag>("njets") );
    pTHToken_ = consumes<float>( HTXSps.getParameter<InputTag>("pTH") );
    pTVToken_ = consumes<float>( HTXSps.getParameter<InputTag>("pTV") );
    newHTXSToken_ = consumes<HTXS::HiggsClassification>( HTXSps.getParameter<InputTag>("ClassificationObj") );

}

DiPhotonJetsTreeMaker::~DiPhotonJetsTreeMaker()
{
}

void
DiPhotonJetsTreeMaker::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup )
{

    // ********************************************************************************

    // access edm objects
    Handle<View<flashgg::DiPhotonCandidate> > diphotons;
    iEvent.getByToken( diphotonToken_, diphotons );

    JetCollectionVector Jets( inputTagJets_.size() );
    for( size_t j = 0; j < inputTagJets_.size(); ++j ) {
        iEvent.getByToken( tokenJets_[j], Jets[j] );
    }

    Handle<View<reco::Vertex> > primaryVertices;
    iEvent.getByToken( vertexToken_, primaryVertices );

    Handle<reco::BeamSpot> recoBeamSpotHandle;
    iEvent.getByToken( beamSpotToken_, recoBeamSpotHandle );

    Handle<double>  rho;
    iEvent.getByToken(rhoTaken_,rho);

    Handle<View<reco::GenParticle> > genParticles;
    iEvent.getByToken( genParticleToken_, genParticles );

    Handle<GenEventInfoProduct> genEventInfo;
    iEvent.getByToken( genEventInfoToken_, genEventInfo );

    Handle<View< PileupSummaryInfo> > pileupInfo;
    iEvent.getByToken( pileUpToken_, pileupInfo );

    Handle<edm::TriggerResults> triggerHandle;
    iEvent.getByToken( triggerToken_, triggerHandle );

    Handle<int> stage0cat, stage1cat, njets;
    Handle<float> pTH, pTV;
    iEvent.getByToken(stage0catToken_, stage0cat);
    iEvent.getByToken(stage1catToken_,stage1cat);
    iEvent.getByToken(njetsToken_,njets);
    iEvent.getByToken(pTHToken_,pTH);
    iEvent.getByToken(pTVToken_,pTV);

    Handle<HTXS::HiggsClassification> htxsClassification;
    iEvent.getByToken(newHTXSToken_,htxsClassification);

    // ********************************************************************************

    initEventStructure();

    if (iEvent.isRealData()) {
        diPhotonJetsInfo.passTrigger = 1;
    } else {
        const edm::TriggerNames &triggerNames = iEvent.triggerNames( *triggerHandle );

        map<string, int> triggerIndices;
        for( unsigned int i = 0; i < triggerNames.triggerNames().size(); i++ ) {
            std::string trimmedName = HLTConfigProvider::removeVersion( triggerNames.triggerName( i ) );
            triggerIndices.emplace(trimmedName, triggerNames.triggerIndex( triggerNames.triggerName( i ) ));
        }

        for (const auto& it : triggerIndices) {
            if (triggerHandle->accept(it.second)) {
                if (it.first == pathName_) {diPhotonJetsInfo.passTrigger = 1; break;}
            }
        }
    }

    diPhotonJetsInfo.rho     = *rho;
    diPhotonJetsInfo.pvz     = primaryVertices->ptrAt(0)->z();
    diPhotonJetsInfo.nvertex = primaryVertices->size();
    if( recoBeamSpotHandle.isValid() ) {
        diPhotonJetsInfo.BSsigmaz = recoBeamSpotHandle->sigmaZ();
    }

    // diphoton loop
    const std::vector<edm::Ptr<flashgg::DiPhotonCandidate> > diphotonsPtrs = diphotons->ptrs();

    if (diphotonsPtrs.size() > 0) {

        Ptr<flashgg::DiPhotonCandidate> diphoPtr = diphotonsPtrs[0];
        diPhotonJetsInfo.dipho_mass                 = diphoPtr->mass();
        diPhotonJetsInfo.dipho_pt                   = diphoPtr->pt();
        diPhotonJetsInfo.dipho_leadPt               = diphoPtr->leadingPhoton()->pt();
        diPhotonJetsInfo.dipho_leadEta              = diphoPtr->leadingPhoton()->eta();
        diPhotonJetsInfo.dipho_leadPhi              = diphoPtr->leadingPhoton()->phi();
        diPhotonJetsInfo.dipho_leadE                = diphoPtr->leadingPhoton()->energy();
        diPhotonJetsInfo.dipho_leadsigEOverE        = diphoPtr->leadingPhoton()->sigEOverE();
        diPhotonJetsInfo.dipho_leadR9               = diphoPtr->leadingPhoton()->full5x5_r9();
        diPhotonJetsInfo.dipho_leadsieie            = diphoPtr->leadingPhoton()->full5x5_sigmaIetaIeta();
        diPhotonJetsInfo.dipho_leadhoe              = diphoPtr->leadingPhoton()->hadronicOverEm();
        diPhotonJetsInfo.dipho_leadIDMVA            = diphoPtr->leadingView()->phoIdMvaWrtChosenVtx();
        diPhotonJetsInfo.dipho_leadGenMatch         = (int) diphoPtr->leadingPhoton()->hasMatchedGenPhoton();
        diPhotonJetsInfo.dipho_leadGenMatchType     = diphoPtr->leadingPhoton()->genMatchType();//enum mcMatch_t { kUnkown = 0, kPrompt, kFake  };
        diPhotonJetsInfo.dipho_subleadPt            = diphoPtr->subLeadingPhoton()->pt();
        diPhotonJetsInfo.dipho_subleadEta           = diphoPtr->subLeadingPhoton()->eta();
        diPhotonJetsInfo.dipho_subleadPhi           = diphoPtr->subLeadingPhoton()->phi();
        diPhotonJetsInfo.dipho_subleadE             = diphoPtr->subLeadingPhoton()->energy();
        diPhotonJetsInfo.dipho_subleadsigEOverE     = diphoPtr->subLeadingPhoton()->sigEOverE();
        diPhotonJetsInfo.dipho_subleadR9            = diphoPtr->subLeadingPhoton()->full5x5_r9();
        diPhotonJetsInfo.dipho_subleadsieie         = diphoPtr->subLeadingPhoton()->full5x5_sigmaIetaIeta();
        diPhotonJetsInfo.dipho_subleadhoe           = diphoPtr->subLeadingPhoton()->hadronicOverEm();
        diPhotonJetsInfo.dipho_subleadIDMVA         = diphoPtr->subLeadingView()->phoIdMvaWrtChosenVtx();
        diPhotonJetsInfo.dipho_subleadGenMatch      = (int) diphoPtr->subLeadingPhoton()->hasMatchedGenPhoton();
        diPhotonJetsInfo.dipho_subleadGenMatchType  = diphoPtr->subLeadingPhoton()->genMatchType();

        diPhotonJetsInfo.selectedvz = diphoPtr->vtx()->position().z();
        diPhotonJetsInfo.genvz = diphoPtr->genPV().z();

        unsigned int jetCollectionIndex = diphoPtr->jetCollectionIndex();
        int Njet = 0;

        for(unsigned int jetLoop = 0; jetLoop < Jets[jetCollectionIndex]->size(); jetLoop++) {
            Ptr<flashgg::Jet> jet = Jets[jetCollectionIndex]->ptrAt(jetLoop);
            if ( !jet->passesJetID(flashgg::Tight2017) ) continue;
            if( fabs( jet->eta() ) > 4.7 ) { continue; }

            float dPhi = deltaPhi( jet->phi(), diPhotonJetsInfo.dipho_leadPhi );
            float dEta = jet->eta() - diPhotonJetsInfo.dipho_leadEta;
            if( sqrt( dPhi * dPhi + dEta * dEta ) < 0.4 ) { continue; }

            dPhi = deltaPhi( jet->phi(), diPhotonJetsInfo.dipho_subleadPhi );
            dEta = jet->eta() - diPhotonJetsInfo.dipho_subleadEta;
            if( sqrt( dPhi * dPhi + dEta * dEta ) < 0.4 ) { continue; }

            Njet++;

            diPhotonJetsInfo.jets_Pt.emplace_back  ( jet->pt()   );
            diPhotonJetsInfo.jets_Eta.emplace_back ( jet->eta()  );
            diPhotonJetsInfo.jets_Phi.emplace_back ( jet->phi()  );
            diPhotonJetsInfo.jets_Mass.emplace_back( jet->mass() );

            diPhotonJetsInfo.jets_QGL.emplace_back( jet->QGL() );
            diPhotonJetsInfo.jets_RMS.emplace_back( jet->rms() );
            diPhotonJetsInfo.jets_puJetIdMVA.emplace_back( jet->puJetIdMVA() );
            diPhotonJetsInfo.jets_GenMatch.emplace_back( (int)jet->hasGenMatch() );

        }//jet loop
        diPhotonJetsInfo.jets_size = Njet;

    }

    if (!iEvent.isRealData()) {
        for( unsigned int PVI = 0; PVI < pileupInfo->size(); ++PVI ) {
            Int_t pu_bunchcrossing = pileupInfo->ptrAt( PVI )->getBunchCrossing();
            if( pu_bunchcrossing == 0 ) {
                diPhotonJetsInfo.nPu = pileupInfo->ptrAt( PVI )->getTrueNumInteractions();
                break;
            }
        }

        if(doHTXS_) {
            diPhotonJetsInfo.genweight     = genEventInfo->weight();
            diPhotonJetsInfo.HTXSstage0cat = htxsClassification->stage0_cat;
            diPhotonJetsInfo.HTXSstage1cat = htxsClassification->stage1_cat_pTjet30GeV;
            diPhotonJetsInfo.HTXSnjets     = htxsClassification->jets30.size();
            diPhotonJetsInfo.HTXSpTH       = htxsClassification->p4decay_higgs.pt();
            diPhotonJetsInfo.HTXSpTV       = htxsClassification->p4decay_V.pt();
        }
        //const std::vector<edm::Ptr<reco::GenParticle> > genParticlesPtrs = genParticles->ptrs();
    }

    DiPhotonJetsTree->Fill();

}


void
DiPhotonJetsTreeMaker::beginJob()
{
    DiPhotonJetsTree = fs_->make<TTree>( "DiPhotonJetsTree", "DiPhotonJets tree" );

    DiPhotonJetsTree->Branch( "nPu"              , &diPhotonJetsInfo.nPu             , "nPu/I"            );
    DiPhotonJetsTree->Branch( "genvz"            , &diPhotonJetsInfo.genvz           , "genvz/F"          );
    DiPhotonJetsTree->Branch( "genweight"        , &diPhotonJetsInfo.genweight       , "genweight/F"      );

    DiPhotonJetsTree->Branch( "rho"              , &diPhotonJetsInfo.rho             , "rho/F"            );
    DiPhotonJetsTree->Branch( "pvz"              , &diPhotonJetsInfo.pvz             , "pvz/F"            );
    DiPhotonJetsTree->Branch( "selectedvz"       , &diPhotonJetsInfo.selectedvz      , "selectedvz/F"     );
    DiPhotonJetsTree->Branch( "nvertex"          , &diPhotonJetsInfo.nvertex         , "nvertex/I"        );
    DiPhotonJetsTree->Branch( "BSsigmaz"         , &diPhotonJetsInfo.BSsigmaz        , "BSsigmaz/F"       );
    DiPhotonJetsTree->Branch( "passTrigger"      , &diPhotonJetsInfo.passTrigger     , "passTrigger/I"    );

    DiPhotonJetsTree->Branch( "dipho_mass"                    , &diPhotonJetsInfo.dipho_mass                   , "dipho_mass/F"                 );
    DiPhotonJetsTree->Branch( "dipho_pt"                      , &diPhotonJetsInfo.dipho_pt                     , "dipho_pt/F"                   );
    DiPhotonJetsTree->Branch( "dipho_leadPt"                  , &diPhotonJetsInfo.dipho_leadPt                 , "dipho_leadPt/F"               );
    DiPhotonJetsTree->Branch( "dipho_leadEta"                 , &diPhotonJetsInfo.dipho_leadEta                , "dipho_leadEta/F"              );
    DiPhotonJetsTree->Branch( "dipho_leadPhi"                 , &diPhotonJetsInfo.dipho_leadPhi                , "dipho_leadPhi/F"              );
    DiPhotonJetsTree->Branch( "dipho_leadE"                   , &diPhotonJetsInfo.dipho_leadE                  , "dipho_leadE/F"                );
    DiPhotonJetsTree->Branch( "dipho_leadsigEOverE"           , &diPhotonJetsInfo.dipho_leadsigEOverE          , "dipho_leadsigEOverE/F"        );
    DiPhotonJetsTree->Branch( "dipho_leadR9"                  , &diPhotonJetsInfo.dipho_leadR9                 , "dipho_leadR9/F"               );
    DiPhotonJetsTree->Branch( "dipho_leadsieie"               , &diPhotonJetsInfo.dipho_leadsieie              , "dipho_leadsieie/F"            );
    DiPhotonJetsTree->Branch( "dipho_leadhoe"                 , &diPhotonJetsInfo.dipho_leadhoe                , "dipho_leadhoe/F"              );
    DiPhotonJetsTree->Branch( "dipho_leadIDMVA"               , &diPhotonJetsInfo.dipho_leadIDMVA              , "dipho_leadIDMVA/F"            );
    DiPhotonJetsTree->Branch( "dipho_leadGenMatch"            , &diPhotonJetsInfo.dipho_leadGenMatch           , "dipho_leadGenMatch/I"         );
    DiPhotonJetsTree->Branch( "dipho_leadGenMatchType"        , &diPhotonJetsInfo.dipho_leadGenMatchType       , "dipho_leadGenMatchType/I"     );
    DiPhotonJetsTree->Branch( "dipho_subleadPt"               , &diPhotonJetsInfo.dipho_subleadPt              , "dipho_subleadPt/F"            );
    DiPhotonJetsTree->Branch( "dipho_subleadEta"              , &diPhotonJetsInfo.dipho_subleadEta             , "dipho_subleadEta/F"           );
    DiPhotonJetsTree->Branch( "dipho_subleadPhi"              , &diPhotonJetsInfo.dipho_subleadPhi             , "dipho_subleadPhi/F"           );
    DiPhotonJetsTree->Branch( "dipho_subleadE"                , &diPhotonJetsInfo.dipho_subleadE               , "dipho_subleadE/F"             );
    DiPhotonJetsTree->Branch( "dipho_subleadsigEOverE"        , &diPhotonJetsInfo.dipho_subleadsigEOverE       , "dipho_subleadsigEOverE/F"     );
    DiPhotonJetsTree->Branch( "dipho_subleadR9"               , &diPhotonJetsInfo.dipho_subleadR9              , "dipho_subleadR9/F"            );
    DiPhotonJetsTree->Branch( "dipho_subleadsieie"            , &diPhotonJetsInfo.dipho_subleadsieie           , "dipho_subleadsieie/F"         );
    DiPhotonJetsTree->Branch( "dipho_subleadhoe"              , &diPhotonJetsInfo.dipho_subleadhoe             , "dipho_subleadhoe/F"           );
    DiPhotonJetsTree->Branch( "dipho_subleadIDMVA"            , &diPhotonJetsInfo.dipho_subleadIDMVA           , "dipho_subleadIDMVA/F"         );
    DiPhotonJetsTree->Branch( "dipho_subleadGenMatch"         , &diPhotonJetsInfo.dipho_subleadGenMatch        , "dipho_subleadGenMatch/I"      );
    DiPhotonJetsTree->Branch( "dipho_subleadGenMatchType"     , &diPhotonJetsInfo.dipho_subleadGenMatchType    , "dipho_subleadGenMatchType/I"  );

    DiPhotonJetsTree->Branch( "jets_size"        , &diPhotonJetsInfo.jets_size      , "jets_size/I" );
    DiPhotonJetsTree->Branch( "jets_Pt"          , &diPhotonJetsInfo.jets_Pt        );
    DiPhotonJetsTree->Branch( "jets_Eta"         , &diPhotonJetsInfo.jets_Eta       );
    DiPhotonJetsTree->Branch( "jets_Phi"         , &diPhotonJetsInfo.jets_Phi       );
    DiPhotonJetsTree->Branch( "jets_Mass"        , &diPhotonJetsInfo.jets_Mass      );
    DiPhotonJetsTree->Branch( "jets_QGL"         , &diPhotonJetsInfo.jets_QGL       );
    DiPhotonJetsTree->Branch( "jets_RMS"         , &diPhotonJetsInfo.jets_RMS       );
    DiPhotonJetsTree->Branch( "jets_puJetIdMVA"  , &diPhotonJetsInfo.jets_puJetIdMVA);
    DiPhotonJetsTree->Branch( "jets_GenMatch"    , &diPhotonJetsInfo.jets_GenMatch  );

    DiPhotonJetsTree->Branch( "HTXSstage0cat"   , &diPhotonJetsInfo.HTXSstage0cat  , "HTXSstage0cat/I" ); 
    DiPhotonJetsTree->Branch( "HTXSstage1cat"   , &diPhotonJetsInfo.HTXSstage1cat  , "HTXSstage1cat/I" ); 
    DiPhotonJetsTree->Branch( "HTXSnjets"       , &diPhotonJetsInfo.HTXSnjets      , "HTXSnjets/I"     ); 
    DiPhotonJetsTree->Branch( "HTXSpTH"         , &diPhotonJetsInfo.HTXSpTH        , "HTXSpTH/F"       ); 
    DiPhotonJetsTree->Branch( "HTXSpTV"         , &diPhotonJetsInfo.HTXSpTV        , "HTXSpTV/F"       ); 

}

void
DiPhotonJetsTreeMaker::endJob()
{
}

void
DiPhotonJetsTreeMaker::initEventStructure()
{

    diPhotonJetsInfo.nPu             = -999;
    diPhotonJetsInfo.genvz           = -999.;
    diPhotonJetsInfo.genweight       = -999.;

    diPhotonJetsInfo.pvz             = -999.;
    diPhotonJetsInfo.selectedvz      = -999.;
    diPhotonJetsInfo.nvertex         = -999;
    diPhotonJetsInfo.BSsigmaz        = -999.;
    diPhotonJetsInfo.passTrigger     = 0;

    diPhotonJetsInfo.dipho_mass                     = -999.;
    diPhotonJetsInfo.dipho_pt                       = -999.;
    diPhotonJetsInfo.dipho_leadPt                   = -999.; 
    diPhotonJetsInfo.dipho_leadEta                  = -999.; 
    diPhotonJetsInfo.dipho_leadPhi                  = -999.; 
    diPhotonJetsInfo.dipho_leadE                    = -999.; 
    diPhotonJetsInfo.dipho_leadsigEOverE            = -999.; 
    diPhotonJetsInfo.dipho_leadR9                   = -999.; 
    diPhotonJetsInfo.dipho_leadsieie                = -999.; 
    diPhotonJetsInfo.dipho_leadhoe                  = -999.; 
    diPhotonJetsInfo.dipho_leadIDMVA                = -999.; 
    diPhotonJetsInfo.dipho_leadGenMatch             = -999; 
    diPhotonJetsInfo.dipho_leadGenMatchType         = -999; 
    diPhotonJetsInfo.dipho_subleadPt                = -999.; 
    diPhotonJetsInfo.dipho_subleadEta               = -999.; 
    diPhotonJetsInfo.dipho_subleadPhi               = -999.; 
    diPhotonJetsInfo.dipho_subleadE                 = -999.; 
    diPhotonJetsInfo.dipho_subleadsigEOverE         = -999.; 
    diPhotonJetsInfo.dipho_subleadR9                = -999.; 
    diPhotonJetsInfo.dipho_subleadsieie             = -999.; 
    diPhotonJetsInfo.dipho_subleadhoe               = -999.; 
    diPhotonJetsInfo.dipho_subleadIDMVA             = -999.; 
    diPhotonJetsInfo.dipho_subleadGenMatch          = -999; 
    diPhotonJetsInfo.dipho_subleadGenMatchType      = -999; 
    
    diPhotonJetsInfo.jets_size   = -999;
    diPhotonJetsInfo.jets_Pt        .clear(); 
    diPhotonJetsInfo.jets_Eta       .clear(); 
    diPhotonJetsInfo.jets_Phi       .clear(); 
    diPhotonJetsInfo.jets_Mass      .clear(); 
    diPhotonJetsInfo.jets_QGL       .clear(); 
    diPhotonJetsInfo.jets_RMS       .clear(); 
    diPhotonJetsInfo.jets_puJetIdMVA.clear(); 
    diPhotonJetsInfo.jets_GenMatch  .clear();

    diPhotonJetsInfo.HTXSstage0cat = -999; 
    diPhotonJetsInfo.HTXSstage1cat = -999; 
    diPhotonJetsInfo.HTXSnjets     = -999; 
    diPhotonJetsInfo.HTXSpTH       = -999.; 
    diPhotonJetsInfo.HTXSpTV       = -999.; 

}

void
DiPhotonJetsTreeMaker::fillDescriptions( edm::ConfigurationDescriptions &descriptions )
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault( desc );
}


DEFINE_FWK_MODULE( DiPhotonJetsTreeMaker );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
