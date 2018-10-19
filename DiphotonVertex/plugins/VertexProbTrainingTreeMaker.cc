
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
#include "DataFormats/Common/interface/PtrVector.h"

#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TTree.h"

// **********************************************************************

// define the structures used to create tree branches and fill the trees

struct GlobalInfo {
    int nPu; 
};

struct DiPhoInfo {

    int ndipho               ;
    int dipho_index          ;
    float dipho_mass         ;
    float pt                 ;
    float dipho_PtLead       ; 
    float dipho_EtaLead      ; 
    float dipho_EtaLeadSC    ; 
    float dipho_R9Lead       ; 
    float dipho_PtSubLead    ; 
    float dipho_EtaSubLead   ; 
    float dipho_EtaSubLeadSC ; 
    float dipho_R9SubLead    ; 

    float LogSumPt2    ;
    float PtBal        ;
    float PtAsym       ;
    float NConv        ;
    float PullConv     ;
    float MVA0         ;
    float MVA1         ;
    float MVA2         ;
    float DZ1          ;
    float DZ2          ;
    float VtxProb      ;

    int nPu            ; 
    float nvertex      ;
    float genvz        ;
    float DZtrue       ;
    float BSsigmaz     ;
    float genWeight    ;

};


// **********************************************************************

using namespace std;
using namespace edm;
using namespace flashgg;


// **********************************************************************

class VertexProbTrainingTreeMaker : public edm::EDAnalyzer
{
public:
    explicit VertexProbTrainingTreeMaker( const edm::ParameterSet & );
    ~VertexProbTrainingTreeMaker();

    static void fillDescriptions( edm::ConfigurationDescriptions &descriptions );


private:

    edm::Service<TFileService> fs_;
    virtual void beginJob() override;
    virtual void analyze( const edm::Event &, const edm::EventSetup & ) override;
    virtual void endJob() override;

    void initEventStructure();
    int getMCTruthVertexIndex( const PtrVector<reco::GenParticle> &gens, const edm::PtrVector<reco::Vertex> &vertices, const double &dzMatch = 0.1 );
    int getSortedIndex( const unsigned int &trueVtxIndex, const unsigned int &sizemax, const Ptr<flashgg::DiPhotonCandidate> &diphoPtr );

    edm::EDGetTokenT<edm::View<flashgg::DiPhotonCandidate> > diphotonToken_;
    edm::EDGetTokenT<edm::View<reco::Vertex> >               vertexToken_;

    edm::EDGetTokenT<reco::BeamSpot > beamSpotToken_;
    edm::EDGetTokenT<double> rhoTaken_;
    edm::EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
    edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;
    edm::EDGetTokenT<View<PileupSummaryInfo>> PileUpToken_;

    bool isZeroVertex_;

    TTree *GlobalTree;
    TTree *DiphoTree;

    GlobalInfo globalInfo;
    DiPhoInfo diphoInfo;

};

// ******************************************************************************************


//
// constructors and destructor
//
VertexProbTrainingTreeMaker::VertexProbTrainingTreeMaker( const edm::ParameterSet &iConfig ):
    diphotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
    vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
    beamSpotToken_( consumes<reco::BeamSpot >( iConfig.getParameter<InputTag>( "BeamSpotTag" ) ) ),
    rhoTaken_( consumes<double>( iConfig.getParameter<InputTag>( "rhoTag" ) ) ),
    genParticleToken_( consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag> ( "GenParticleTag" ) ) ),
    genEventInfoToken_( consumes<GenEventInfoProduct>( iConfig.getParameter <edm::InputTag> ( "GenEventInfoTag" ) ) ),
    PileUpToken_( consumes<View<PileupSummaryInfo> >( iConfig.getParameter<InputTag> ( "PileUpTag" ) ) )
{
    isZeroVertex_ = iConfig.getParameter<bool>( "isZeroVertex" );
}


VertexProbTrainingTreeMaker::~VertexProbTrainingTreeMaker()
{
}

void
VertexProbTrainingTreeMaker::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup )
{

    // ********************************************************************************

    // access edm objects
    Handle<View<flashgg::DiPhotonCandidate> > diphotons;
    iEvent.getByToken( diphotonToken_, diphotons );
    const std::vector<edm::Ptr<flashgg::DiPhotonCandidate> > &diphotonPointers = diphotons->ptrs();

    Handle<View<reco::Vertex> > primaryVertices;
    iEvent.getByToken( vertexToken_, primaryVertices );
    const std::vector<edm::Ptr<reco::Vertex> > &pvPointers = primaryVertices->ptrs();

    Handle<reco::BeamSpot> recoBeamSpotHandle;
    iEvent.getByToken( beamSpotToken_, recoBeamSpotHandle );

    Handle<double>  rho;
    iEvent.getByToken(rhoTaken_,rho);

    Handle<View<reco::GenParticle> > genParticles;
    iEvent.getByToken( genParticleToken_, genParticles );
    //const std::vector<edm::Ptr<reco::GenParticle> > &gens = genParticles->ptrs();
    
    Handle<GenEventInfoProduct> genEventInfo;
    iEvent.getByToken( genEventInfoToken_, genEventInfo );

    Handle<View< PileupSummaryInfo> > PileupInfo;
    iEvent.getByToken( PileUpToken_, PileupInfo ); 


    // ********************************************************************************

    initEventStructure();

    // pileup info
    for( unsigned int PVI = 0; PVI < PileupInfo->size(); ++PVI ) {
        Int_t pu_bunchcrossing = PileupInfo->ptrAt( PVI )->getBunchCrossing();
        if( pu_bunchcrossing == 0 ) {
            globalInfo.nPu = PileupInfo->ptrAt( PVI )->getTrueNumInteractions();
            break;
        }
    }
    GlobalTree->Fill();

    double BSsigmaz = 0;
    if( recoBeamSpotHandle.isValid() ) {
        BSsigmaz = recoBeamSpotHandle->sigmaZ();
    }

    for( size_t idipho = 0; idipho < diphotonPointers.size(); idipho++ ) {

        Ptr<flashgg::DiPhotonCandidate> diphoPtr = diphotonPointers[idipho];

        diphoInfo.ndipho = diphotonPointers.size();
        diphoInfo.dipho_index = idipho;
        diphoInfo.dipho_mass = diphoPtr->mass();
        diphoInfo.pt = diphoPtr->pt();

        diphoInfo.dipho_PtLead       = diphoPtr->leadingPhoton()->pt();
        diphoInfo.dipho_EtaLead      = diphoPtr->leadingPhoton()->eta();
        diphoInfo.dipho_EtaLeadSC    = diphoPtr->leadingPhoton()->superCluster()->eta();
        diphoInfo.dipho_R9Lead       = diphoPtr->leadingPhoton()->full5x5_r9();
        diphoInfo.dipho_PtSubLead    = diphoPtr->subLeadingPhoton()->pt();
        diphoInfo.dipho_EtaSubLead   = diphoPtr->subLeadingPhoton()->eta();
        diphoInfo.dipho_EtaSubLeadSC = diphoPtr->subLeadingPhoton()->superCluster()->eta();
        diphoInfo.dipho_R9SubLead    = diphoPtr->subLeadingPhoton()->full5x5_r9();

        if (!isZeroVertex_) {
            diphoInfo.LogSumPt2 = diphoPtr->logSumPt2();
            diphoInfo.PtBal  =  diphoPtr->ptBal();
            diphoInfo.PtAsym  =  diphoPtr->ptAsym();
            diphoInfo.NConv  =  diphoPtr->nConv();
            diphoInfo.PullConv  =  diphoPtr->pullConv();
            diphoInfo.MVA0 = diphoPtr->mva0();
            diphoInfo.MVA1 = diphoPtr->mva1();
            diphoInfo.MVA2 = diphoPtr->mva2();
            diphoInfo.DZ1 = diphoPtr->dZ1();
            diphoInfo.DZ2 = diphoPtr->dZ2();
            diphoInfo.VtxProb = diphoPtr->vtxProbMVA();
        }

        diphoInfo.nPu = globalInfo.nPu;
        diphoInfo.nvertex = pvPointers.size();
        diphoInfo.genvz = diphoPtr->genPV().z();
        diphoInfo.DZtrue = fabs( diphoPtr->vtx()->position().z() - diphoInfo.genvz );
        diphoInfo.BSsigmaz = BSsigmaz;
        diphoInfo.genWeight = genEventInfo->weight() > 0. ? 1. : -1.;


        DiphoTree->Fill();

    }  // end diphoton candidate loop

}

void
VertexProbTrainingTreeMaker::beginJob()
{
    GlobalTree = fs_->make<TTree>( "GlobalTree", "" );
    GlobalTree->Branch( "nPu"                 , &globalInfo.nPu               , "nPu/I"                );

    DiphoTree = fs_->make<TTree>( "DiphoTree", "per-diphoton tree" );
    DiphoTree->Branch( "nPu"                  , &diphoInfo.nPu                , "nPu/I"                );
    DiphoTree->Branch( "ndipho"               , &diphoInfo.ndipho             , "ndipho/I"             );
    DiphoTree->Branch( "dipho_index"          , &diphoInfo.dipho_index        , "dipho_index/I"        );
    DiphoTree->Branch( "dipho_mass"           , &diphoInfo.dipho_mass         , "dipho_mass/F"         );
    DiphoTree->Branch( "pt"                   , &diphoInfo.pt                 , "pt/F"                 );
    DiphoTree->Branch( "dipho_PtLead"         , &diphoInfo.dipho_PtLead       , "dipho_PtLead/F"       );
    DiphoTree->Branch( "dipho_EtaLead"        , &diphoInfo.dipho_EtaLead      , "dipho_EtaLead/F"      );
    DiphoTree->Branch( "dipho_EtaLeadSC"      , &diphoInfo.dipho_EtaLeadSC    , "dipho_EtaLeadSC/F"    );
    DiphoTree->Branch( "dipho_R9Lead"         , &diphoInfo.dipho_R9Lead       , "dipho_R9Lead/F"       );
    DiphoTree->Branch( "dipho_PtSubLead"      , &diphoInfo.dipho_PtSubLead    , "dipho_PtSubLead/F"    );
    DiphoTree->Branch( "dipho_EtaSubLead"     , &diphoInfo.dipho_EtaSubLead   , "dipho_EtaSubLead/F"   );
    DiphoTree->Branch( "dipho_EtaSubLeadSC"   , &diphoInfo.dipho_EtaSubLeadSC , "dipho_EtaSubLeadSC/F" );
    DiphoTree->Branch( "dipho_R9SubLead"      , &diphoInfo.dipho_R9SubLead    , "dipho_R9SubLead/F"    );
    DiphoTree->Branch( "LogSumPt2"            , &diphoInfo.LogSumPt2          , "LogSumPt2/F"          );
    DiphoTree->Branch( "PtBal"                , &diphoInfo.PtBal              , "PtBal/F"              );
    DiphoTree->Branch( "PtAsym"               , &diphoInfo.PtAsym             , "PtAsym/F"             );
    DiphoTree->Branch( "NConv"                , &diphoInfo.NConv              , "NConv/F"              );
    DiphoTree->Branch( "PullConv"             , &diphoInfo.PullConv           , "PullConv/F"           );
    DiphoTree->Branch( "MVA0"                 , &diphoInfo.MVA0               , "MVA0/F"               );
    DiphoTree->Branch( "MVA1"                 , &diphoInfo.MVA1               , "MVA1/F"               );
    DiphoTree->Branch( "MVA2"                 , &diphoInfo.MVA2               , "MVA2/F"               );
    DiphoTree->Branch( "DZ1"                  , &diphoInfo.DZ1                , "DZ1/F"                );
    DiphoTree->Branch( "DZ2"                  , &diphoInfo.DZ2                , "DZ2/F"                );
    DiphoTree->Branch( "vtxProb"              , &diphoInfo.VtxProb            , "vtxProb/F"            );
    DiphoTree->Branch( "nPu"                  , &diphoInfo.nPu                , "nPu/I"                );
    DiphoTree->Branch( "NVert"                , &diphoInfo.nvertex            , "NVert/F"              );
    DiphoTree->Branch( "DZtrue"               , &diphoInfo.DZtrue             , "DZtrue/F"             );
    DiphoTree->Branch( "BSsigmaz"             , &diphoInfo.BSsigmaz           , "BSsigmaz/F"           );
    DiphoTree->Branch( "genWeight"            , &diphoInfo.genWeight          , "genWeight/F"          );
}

void
VertexProbTrainingTreeMaker::endJob()
{
}

void
VertexProbTrainingTreeMaker::initEventStructure()
{
    globalInfo.nPu                = -999;

    diphoInfo.nPu                 = -999; 
    diphoInfo.ndipho              = -999; 
    diphoInfo.dipho_index         = -999; 
    diphoInfo.dipho_mass          = -999; 
    diphoInfo.pt                  = -999; 
    diphoInfo.dipho_PtLead        = -999; 
    diphoInfo.dipho_EtaLead       = -999; 
    diphoInfo.dipho_EtaLeadSC     = -999; 
    diphoInfo.dipho_R9Lead        = -999; 
    diphoInfo.dipho_PtSubLead     = -999; 
    diphoInfo.dipho_EtaSubLead    = -999; 
    diphoInfo.dipho_EtaSubLeadSC  = -999; 
    diphoInfo.dipho_R9SubLead     = -999; 
    diphoInfo.LogSumPt2           = -999; 
    diphoInfo.PtBal               = -999; 
    diphoInfo.PtAsym              = -999; 
    diphoInfo.NConv               = -999; 
    diphoInfo.PullConv            = -999; 
    diphoInfo.MVA0                = -999; 
    diphoInfo.MVA1                = -999; 
    diphoInfo.MVA2                = -999; 
    diphoInfo.DZ1                 = -999; 
    diphoInfo.DZ2                 = -999; 
    diphoInfo.VtxProb             = -999; 
    diphoInfo.nPu                 = -999; 
    diphoInfo.nvertex             = -999; 
    diphoInfo.DZtrue              = -999; 
    diphoInfo.BSsigmaz            = -999; 
    diphoInfo.genWeight           = -999; 
}

void
VertexProbTrainingTreeMaker::fillDescriptions( edm::ConfigurationDescriptions &descriptions )
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault( desc );
}

int 
VertexProbTrainingTreeMaker::getSortedIndex( const unsigned int &trueVtxIndex, const unsigned int &sizemax, const Ptr<flashgg::DiPhotonCandidate> &diphoPtr )
{

    for( unsigned int j = 0; j < sizemax; j++ ) {
        int index = diphoPtr->mvaSortedIndex( j );
        if( index < 0 ) { continue; }
        if( ( unsigned int ) index == trueVtxIndex ) { return j; }
    }
    return -1;
}

int 
VertexProbTrainingTreeMaker::getMCTruthVertexIndex( const PtrVector<reco::GenParticle> &gens , const PtrVector<reco::Vertex> &vertices, const double &dzMatch )
{

    reco::Vertex::Point hardVertex( 0, 0, 0 );

    for( unsigned int genLoop = 0 ; genLoop < gens.size(); genLoop++ ) {

        if( fabs( gens[genLoop]->pdgId() ) < 10 || fabs( gens[genLoop]->pdgId() ) == 25 ) {
            hardVertex.SetCoordinates( gens[genLoop]->vx(), gens[genLoop]->vy(), gens[genLoop]->vz() );
            break;
        }
    }

    int  ivMatch = 0;
    double dzMin = 999;
    for( unsigned int iv = 0; iv < vertices.size(); iv++ ) {
        double dz = fabs( vertices[iv]->z() - hardVertex.z() );
        if( dz < dzMin ) {
            ivMatch = iv;
            dzMin   = dz;
        }
    }


    if( dzMin < dzMatch ) { return ivMatch; }

    return -1;
}


DEFINE_FWK_MODULE( VertexProbTrainingTreeMaker );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
