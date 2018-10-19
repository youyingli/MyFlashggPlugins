
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

#include "flashgg/DataFormats/interface/PhotonJetCandidate.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TTree.h"

// **********************************************************************

// define the structures used to create tree branches and fill the trees


struct GlobalInfo {
    
    int nvtx;
    int nPu;
    int nPuTrue;
    float BSsigmaz;
    int nphoj;
    float evWeight;

};

struct PhotonJetInfo {

    int nvtx;
    int nPu;
    int nPuTrue;
    float BSsigmaz;
    float genweight;
    float evWeight;

    int   JetTagVtxIndex; 
    float JetTagVtxIDMva; 
    int   UnJetTagVtxIndex;
    float UnJetTagVtxIDMva;
    int   SelectedvtxIndex;
    float SelectedVtxIDMva;
    int   UnSelectedvtxIndex;
    float UnSelectedVtxIDMva;
    float dZfromRecoPV;    
    float nConv;           
    float vtxProbMVA;      
    float photonjet_Pt;     
    float photon_pt;       
    float photon_eta;       
    float jet_pt;          
    float jet_eta;         
    float photonjet_sumPt;         

    float JetTagVtx_logSumPt2;       
    float JetTagVtx_ptBal;           
    float JetTagVtx_ptAsym;          
    float JetTagVtx_pullConv;       
    float UnJetTagVtx_logSumPt2;       
    float UnJetTagVtx_ptBal;           
    float UnJetTagVtx_ptAsym;          
    float UnJetTagVtx_pullConv;        
    float Selectedvtx_logSumPt2;       
    float Selectedvtx_ptBal;           
    float Selectedvtx_ptAsym;          
    float Selectedvtx_pullConv;       
    float UnSelectedvtx_logSumPt2;       
    float UnSelectedvtx_ptBal;           
    float UnSelectedvtx_ptAsym;          
    float UnSelectedvtx_pullConv;        

};

// **********************************************************************

using namespace std;
using namespace edm;
using namespace flashgg;


// **********************************************************************

class PhotonJetValidationTreeMaker : public edm::EDAnalyzer
{
public:
    explicit PhotonJetValidationTreeMaker( const edm::ParameterSet & );
    ~PhotonJetValidationTreeMaker();

    static void fillDescriptions( edm::ConfigurationDescriptions &descriptions );


private:

    edm::Service<TFileService> fs_;

    virtual void beginJob() override;
    virtual void analyze( const edm::Event &, const edm::EventSetup & ) override;
    virtual void endJob() override;

    void initEventStructure();
    int sortedIndex( const unsigned int &VtxIndex, const unsigned int &sizemax, const Ptr<flashgg::PhotonJetCandidate> &photonjet );

    edm::EDGetTokenT<edm::View<flashgg::PhotonJetCandidate> >       photonJetToken_;
    edm::EDGetTokenT<edm::View<reco::Vertex>>                       vertexToken_;
    
    edm::EDGetTokenT<GenEventInfoProduct>                           genEventInfoToken_;
    edm::EDGetTokenT<reco::BeamSpot>                                beamSpotToken_;
    edm::EDGetTokenT<View<PileupSummaryInfo>>                       PileUpToken_;
    double evWeight_;

    TTree *GlobalTree;
    TTree *PhotonJetTree;
    GlobalInfo globalInfo;
    PhotonJetInfo phojInfo;
};

// ******************************************************************************************


//
// constructors and destructor
//
PhotonJetValidationTreeMaker::PhotonJetValidationTreeMaker( const edm::ParameterSet &iConfig ):
    photonJetToken_( consumes<View<flashgg::PhotonJetCandidate> >( iConfig.getParameter<InputTag> ( "PhotonJetTag" ) ) ),
    vertexToken_( consumes<View<reco::Vertex>>( iConfig.getParameter<InputTag> ( "vertexTag" ) ) ),
    genEventInfoToken_( consumes<GenEventInfoProduct>( iConfig.getParameter<InputTag> ( "GenEventInfo" ) ) ),
    beamSpotToken_( consumes<reco::BeamSpot>( iConfig.getParameter<InputTag>( "BeamSpotTag" ) ) ),
    PileUpToken_( consumes<View<PileupSummaryInfo>>( iConfig.getParameter<InputTag> ( "PileUpTag" ) ) )
{
    evWeight_   = iConfig.getParameter<double>( "evWeight" );
}


PhotonJetValidationTreeMaker::~PhotonJetValidationTreeMaker()
{
}

void
PhotonJetValidationTreeMaker::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup )
{

    // ********************************************************************************

    // access edm objects
    
    Handle<View<flashgg::PhotonJetCandidate> > photonJets;
    iEvent.getByToken(photonJetToken_, photonJets);

    Handle<View<reco::Vertex>> VertexHandle_;
    iEvent.getByToken(vertexToken_, VertexHandle_);
 
    Handle<View<PileupSummaryInfo>> PileUpHandle_;
    iEvent.getByToken(PileUpToken_, PileUpHandle_);

    Handle<GenEventInfoProduct> GenInfoHandle_;
    iEvent.getByToken(genEventInfoToken_, GenInfoHandle_);

    Handle<reco::BeamSpot> recoBeamSpotHandle;
    iEvent.getByToken( beamSpotToken_, recoBeamSpotHandle );
    double BSsigmaz = 0.;
    if( recoBeamSpotHandle.isValid() ) BSsigmaz = recoBeamSpotHandle->sigmaZ();

    // ********************************************************************************

    initEventStructure();

    globalInfo.BSsigmaz = BSsigmaz;
    globalInfo.nvtx  = VertexHandle_->size();
    globalInfo.nphoj = photonJets->size();
    globalInfo.evWeight = evWeight_;

    if (!iEvent.isRealData()) {
        for( unsigned int PVI = 0; PVI < PileUpHandle_->size(); ++PVI ) {
            Int_t pu_bunchcrossing = PileUpHandle_->ptrAt( PVI )->getBunchCrossing();
            if( pu_bunchcrossing == 0 ) {
                globalInfo.nPu = PileUpHandle_->ptrAt( PVI )->getPU_NumInteractions();
                globalInfo.nPuTrue = PileUpHandle_->ptrAt( PVI )->getTrueNumInteractions();
                break;
            }
        }
    }
    GlobalTree->Fill();

    // photon jet loop
    for( size_t iphoj = 0; iphoj < photonJets->size(); iphoj++ ) {

        Ptr<flashgg::PhotonJetCandidate> photonJet = photonJets->ptrAt(iphoj);

        if (photonJet->nConv() != 1) continue;

        phojInfo.BSsigmaz = globalInfo.BSsigmaz;
        phojInfo.nvtx     = globalInfo.nvtx;
        if (!iEvent.isRealData()) {
            phojInfo.nPu      = globalInfo.nPu;
            phojInfo.nPuTrue  = globalInfo.nPuTrue;
            phojInfo.genweight= GenInfoHandle_ -> weight();
            phojInfo.evWeight = evWeight_;
        }

        vector<int> pvVecNoJetTagged;
        for( int i = 0 ; i < (int)VertexHandle_->size() ; i++ ) {
            if( i != photonJet->vertexIndexJet() ) pvVecNoJetTagged.push_back( i ); 
        }

        vector<int> pvVecNoSelected;
        for( int i = 0 ; i < (int)VertexHandle_->size() ; i++ ) {
            if( i != photonJet->vertexIndex() ) pvVecNoSelected.push_back( i ); 
        }

        int irand = -1;
        if( pvVecNoJetTagged.size() > 1 ) { irand = rand() % pvVecNoJetTagged.size(); }
 
        int randNoJetTaggedVtxIndex = -1;
        if( irand != -1 ) { randNoJetTaggedVtxIndex = pvVecNoJetTagged[irand]; }

        int irands = -1;
        if( pvVecNoSelected.size() > 1 ) { irands = rand() % pvVecNoSelected.size(); }

        int randNoSelectedVtxIndex = -1;
        if( irands != -1 ) { randNoSelectedVtxIndex = pvVecNoSelected[irands]; }


        int JetTagVtxIndex = sortedIndex(photonJet->vertexIndexJet(), VertexHandle_->size(), photonJet);
        int UnJetTagVtxIndex = sortedIndex(randNoJetTaggedVtxIndex, VertexHandle_->size(), photonJet);
        int SelectedVtxIndex = sortedIndex(photonJet->vertexIndex(), VertexHandle_->size(), photonJet);
        int UnSelectedVtxIndex = sortedIndex(randNoSelectedVtxIndex, VertexHandle_->size(), photonJet);


        phojInfo.JetTagVtxIndex     = photonJet->vertexIndexJet();
        phojInfo.JetTagVtxIDMva     = photonJet->mva(JetTagVtxIndex);
        phojInfo.UnJetTagVtxIndex   = photonJet->mvaSortedIndex(UnJetTagVtxIndex);
        phojInfo.UnJetTagVtxIDMva   = photonJet->mva(UnJetTagVtxIndex);
        phojInfo.SelectedvtxIndex   = photonJet->vertexIndex();
        phojInfo.SelectedVtxIDMva   = photonJet->mva(SelectedVtxIndex);
        phojInfo.UnSelectedvtxIndex = photonJet->mvaSortedIndex(UnSelectedVtxIndex);
        phojInfo.UnSelectedVtxIDMva = photonJet->mva(UnSelectedVtxIndex);

        phojInfo.dZfromRecoPV     = photonJet->dZfromRecoPV();
        phojInfo.nConv            = photonJet->nConv();
        phojInfo.vtxProbMVA       = photonJet->vtxProbMVA();
        phojInfo.photonjet_Pt     = photonJet->photonjetPt();
        phojInfo.photon_pt        = photonJet->photon()->pt();
        phojInfo.photon_eta       = photonJet->photon()->eta();
        phojInfo.jet_pt           = photonJet->jet()->pt();
        phojInfo.jet_eta          = photonJet->jet()->eta();
        phojInfo.photonjet_sumPt  = photonJet->sumPt();

        phojInfo.JetTagVtx_logSumPt2          = photonJet->logSumPt2(JetTagVtxIndex);
        phojInfo.JetTagVtx_ptBal              = photonJet->ptBal(JetTagVtxIndex);
        phojInfo.JetTagVtx_ptAsym             = photonJet->ptAsym(JetTagVtxIndex);
        phojInfo.JetTagVtx_pullConv           = photonJet->pullConv(JetTagVtxIndex);
        phojInfo.UnJetTagVtx_logSumPt2        = photonJet->logSumPt2(UnJetTagVtxIndex);
        phojInfo.UnJetTagVtx_ptBal            = photonJet->ptBal(UnJetTagVtxIndex);
        phojInfo.UnJetTagVtx_ptAsym           = photonJet->ptAsym(UnJetTagVtxIndex);
        phojInfo.UnJetTagVtx_pullConv         = photonJet->pullConv(UnJetTagVtxIndex);
        phojInfo.Selectedvtx_logSumPt2        = photonJet->logSumPt2(SelectedVtxIndex);
        phojInfo.Selectedvtx_ptBal            = photonJet->ptBal(SelectedVtxIndex);
        phojInfo.Selectedvtx_ptAsym           = photonJet->ptAsym(SelectedVtxIndex);
        phojInfo.Selectedvtx_pullConv         = photonJet->pullConv(SelectedVtxIndex);
        phojInfo.UnSelectedvtx_logSumPt2      = photonJet->logSumPt2(UnSelectedVtxIndex);
        phojInfo.UnSelectedvtx_ptBal          = photonJet->ptBal(UnSelectedVtxIndex);
        phojInfo.UnSelectedvtx_ptAsym         = photonJet->ptAsym(UnSelectedVtxIndex);
        phojInfo.UnSelectedvtx_pullConv       = photonJet->pullConv(UnSelectedVtxIndex);

        PhotonJetTree->Fill();

        break; // Only choose the first candidate

    }  // end photon jet candidate loop

}


void
PhotonJetValidationTreeMaker::beginJob()
{
    GlobalTree = fs_->make<TTree>( "GlobalTree", "global tree" );
    GlobalTree->Branch( "BSsigmaz"  , &globalInfo.BSsigmaz,   "BSsigmaz/F" );
    GlobalTree->Branch( "nvtx"      , &globalInfo.nvtx    ,   "nvtx/I"     );
    GlobalTree->Branch( "nPu"       , &globalInfo.nPu     ,   "nPu/I"      );
    GlobalTree->Branch( "nPuTrue"   , &globalInfo.nPuTrue ,   "nPuTrue/I"  );
    GlobalTree->Branch( "nphoj"     , &globalInfo.nphoj   ,   "nphoj/I"    );
    GlobalTree->Branch( "evWeight"     , &globalInfo.evWeight   ,   "evWeight/F"    );

    PhotonJetTree = fs_->make<TTree>( "PhotonJetTree", "photon-jet tree" );
    PhotonJetTree->Branch( "BSsigmaz"           , &phojInfo.BSsigmaz        ,   "BSsigmaz/F"            );
    PhotonJetTree->Branch( "nvtx"               , &phojInfo.nvtx            ,   "nvtx/I"                );
    PhotonJetTree->Branch( "nPu"                , &phojInfo.nPu             ,   "nPu/I"                 );
    PhotonJetTree->Branch( "nPuTrue"            , &phojInfo.nPuTrue         ,   "nPuTrue/I"             );
    PhotonJetTree->Branch( "genweight"          , &phojInfo.genweight       ,   "genweight/F"           );
    PhotonJetTree->Branch( "evWeight"           , &phojInfo.evWeight        ,   "evWeight/F"            );
    PhotonJetTree->Branch( "JetTagVtxIndex"     , &phojInfo.JetTagVtxIndex  ,   "JetTagVtxIndex/I"      );
    PhotonJetTree->Branch( "JetTagVtxIDMva"     , &phojInfo.JetTagVtxIDMva  ,   "JetTagVtxIDMva/F"      );
    PhotonJetTree->Branch( "UnJetTagVtxIndex"   , &phojInfo.UnJetTagVtxIndex,   "UnJetTagVtxIndex/I"    );
    PhotonJetTree->Branch( "UnJetTagVtxIDMva"   , &phojInfo.UnJetTagVtxIDMva,   "UnJetTagVtxIDMva/F"    );
    PhotonJetTree->Branch( "SelectedvtxIndex"   , &phojInfo.SelectedvtxIndex,   "SelectedvtxIndex/I"    );
    PhotonJetTree->Branch( "SelectedVtxIDMva"   , &phojInfo.SelectedVtxIDMva,   "SelectedVtxIDMva/F"    );
    PhotonJetTree->Branch( "UnSelectedvtxIndex"   , &phojInfo.UnSelectedvtxIndex,   "UnSelectedvtxIndex/I"    );
    PhotonJetTree->Branch( "UnSelectedVtxIDMva"   , &phojInfo.UnSelectedVtxIDMva,   "UnSelectedVtxIDMva/F"    );
    PhotonJetTree->Branch( "dZfromRecoPV"       , &phojInfo.dZfromRecoPV    ,   "dZfromRecoPV/F"        );
    PhotonJetTree->Branch( "nConv"              , &phojInfo.nConv           ,   "nConv/F"               );
    PhotonJetTree->Branch( "vtxProbMVA"         , &phojInfo.vtxProbMVA      ,   "vtxProbMVA/F"          );
    PhotonJetTree->Branch( "photonjet_Pt"       , &phojInfo.photonjet_Pt    ,   "photonjet_Pt/F"        );
    PhotonJetTree->Branch( "photon_pt"          , &phojInfo.photon_pt       ,   "photon_pt/F"           );
    PhotonJetTree->Branch( "photon_eta"         , &phojInfo.photon_eta      ,   "photon_eta/F"          );
    PhotonJetTree->Branch( "jet_pt"             , &phojInfo.jet_pt          ,   "jet_pt/F"              );
    PhotonJetTree->Branch( "jet_eta"            , &phojInfo.jet_eta         ,   "jet_eta/F"             );
    PhotonJetTree->Branch( "photonjet_sumPt"    , &phojInfo.photonjet_sumPt ,   "photonjet_sumPt/F"     );

    PhotonJetTree->Branch( "JetTagVtx_logSumPt2"    , &phojInfo.JetTagVtx_logSumPt2         ,   "JetTagVtx_logSumPt2/F"     );
    PhotonJetTree->Branch( "JetTagVtx_ptBal"        , &phojInfo.JetTagVtx_ptBal             ,   "JetTagVtx_ptBal/F"         );
    PhotonJetTree->Branch( "JetTagVtx_ptAsym"       , &phojInfo.JetTagVtx_ptAsym            ,   "JetTagVtx_ptAsym/F"        );
    PhotonJetTree->Branch( "JetTagVtx_pullConv"     , &phojInfo.JetTagVtx_pullConv          ,   "JetTagVtx_pullConv/F"      );
    PhotonJetTree->Branch( "UnJetTagVtx_logSumPt2"  , &phojInfo.UnJetTagVtx_logSumPt2       ,   "UnJetTagVtx_logSumPt2/F"   );
    PhotonJetTree->Branch( "UnJetTagVtx_ptBal"      , &phojInfo.UnJetTagVtx_ptBal           ,   "UnJetTagVtx_ptBal/F"       );
    PhotonJetTree->Branch( "UnJetTagVtx_ptAsym"     , &phojInfo.UnJetTagVtx_ptAsym          ,   "UnJetTagVtx_ptAsym/F"      );
    PhotonJetTree->Branch( "UnJetTagVtx_pullConv"   , &phojInfo.UnJetTagVtx_pullConv        ,   "UnJetTagVtx_pullConv/F"    );
    PhotonJetTree->Branch( "Selectedvtx_logSumPt2"  , &phojInfo.Selectedvtx_logSumPt2       ,   "Selectedvtx_logSumPt2/F"   );
    PhotonJetTree->Branch( "Selectedvtx_ptBal"      , &phojInfo.Selectedvtx_ptBal           ,   "Selectedvtx_ptBal/F"       );
    PhotonJetTree->Branch( "Selectedvtx_ptAsym"     , &phojInfo.Selectedvtx_ptAsym          ,   "Selectedvtx_ptAsym/F"      );
    PhotonJetTree->Branch( "Selectedvtx_pullConv"   , &phojInfo.Selectedvtx_pullConv        ,   "Selectedvtx_pullConv/F"    );
    PhotonJetTree->Branch( "UnSelectedvtx_logSumPt2"  , &phojInfo.UnSelectedvtx_logSumPt2       ,   "UnSelectedvtx_logSumPt2/F"   );
    PhotonJetTree->Branch( "UnSelectedvtx_ptBal"      , &phojInfo.UnSelectedvtx_ptBal           ,   "UnSelectedvtx_ptBal/F"       );
    PhotonJetTree->Branch( "UnSelectedvtx_ptAsym"     , &phojInfo.UnSelectedvtx_ptAsym          ,   "UnSelectedvtx_ptAsym/F"      );
    PhotonJetTree->Branch( "UnSelectedvtx_pullConv"   , &phojInfo.UnSelectedvtx_pullConv        ,   "UnSelectedvtx_pullConv/F"    );

}

void
PhotonJetValidationTreeMaker::endJob()
{
}

void
PhotonJetValidationTreeMaker::initEventStructure()
{
    globalInfo.BSsigmaz           = -999.; 
    globalInfo.nvtx               = -999.; 
    globalInfo.nPu                = -999.; 
    globalInfo.nPuTrue            = -999.; 
    globalInfo.nphoj              = -999.; 
    globalInfo.evWeight           = 1.; 

    phojInfo.BSsigmaz           = -999.; 
    phojInfo.nvtx               = -999.; 
    phojInfo.nPu                = -999.; 
    phojInfo.nPuTrue            = -999.; 
    phojInfo.genweight          = -999.; 
    phojInfo.evWeight           = 1.; 
    phojInfo.JetTagVtxIndex     = -999.; 
    phojInfo.JetTagVtxIDMva     = -999.; 
    phojInfo.UnJetTagVtxIndex   = -999.;
    phojInfo.UnJetTagVtxIDMva   = -999.;
    phojInfo.SelectedvtxIndex   = -999.;
    phojInfo.SelectedVtxIDMva   = -999.;
    phojInfo.UnSelectedvtxIndex   = -999.;
    phojInfo.UnSelectedVtxIDMva   = -999.;
    phojInfo.dZfromRecoPV       = -999.; 
    phojInfo.nConv              = -999.; 
    phojInfo.vtxProbMVA         = -999.; 
    phojInfo.photonjet_Pt       = -999.; 
    phojInfo.photon_pt          = -999.; 
    phojInfo.photon_eta         = -999.; 
    phojInfo.jet_pt             = -999.;
    phojInfo.jet_eta            = -999.;
    phojInfo.photonjet_sumPt    = -999.;

    phojInfo.JetTagVtx_logSumPt2         = -999.;   
    phojInfo.JetTagVtx_ptBal             = -999.;       
    phojInfo.JetTagVtx_ptAsym            = -999.;      
    phojInfo.JetTagVtx_pullConv          = -999.;    
    phojInfo.UnJetTagVtx_logSumPt2       = -999.; 
    phojInfo.UnJetTagVtx_ptBal           = -999.;     
    phojInfo.UnJetTagVtx_ptAsym          = -999.;    
    phojInfo.UnJetTagVtx_pullConv        = -999.;  
    phojInfo.Selectedvtx_logSumPt2       = -999.; 
    phojInfo.Selectedvtx_ptBal           = -999.;     
    phojInfo.Selectedvtx_ptAsym          = -999.;    
    phojInfo.Selectedvtx_pullConv        = -999.; 
    phojInfo.UnSelectedvtx_logSumPt2       = -999.; 
    phojInfo.UnSelectedvtx_ptBal           = -999.;     
    phojInfo.UnSelectedvtx_ptAsym          = -999.;    
    phojInfo.UnSelectedvtx_pullConv        = -999.;  

}

void
PhotonJetValidationTreeMaker::fillDescriptions( edm::ConfigurationDescriptions &descriptions )
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault( desc );
}

int 
PhotonJetValidationTreeMaker::sortedIndex( const unsigned int &VtxIndex, const unsigned int &sizemax, const Ptr<flashgg::PhotonJetCandidate> &photonjet )
{
    for( unsigned int j = 0; j < sizemax; j++ ) {
        int index = photonjet->mvaSortedIndex( j );
        if( index < 0 ) { continue; }
        if( ( unsigned int ) index == VtxIndex ) { return j; }
    }
    return -1;
}

DEFINE_FWK_MODULE( PhotonJetValidationTreeMaker );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

