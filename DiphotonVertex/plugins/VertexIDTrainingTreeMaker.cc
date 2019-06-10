
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

#include "flashgg/DataFormats/interface/VertexCandidateMap.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "TTree.h"

// **********************************************************************

// define the structures used to create tree branches and fill the trees

struct SignalInfo {

    int nvertex;
    int ndipho;
    int dipho_index;

    float LogSumPt2;
    float PtBal;
    float PtAsym;
    float NConv;
    float PullConv;

    float dipho_mass;
    float genweight;
};

struct BackgroundInfo {

    int nvertex;
    int ndipho;
    int dipho_index;

    float LogSumPt2;
    float PtBal;
    float PtAsym;
    float NConv;
    float PullConv;

    float dipho_mass;
    float genweight;
};

// **********************************************************************

using namespace std;
using namespace edm;
using namespace flashgg;


// **********************************************************************

class VertexIDTrainingTreeMaker : public edm::EDAnalyzer
{
public:
    explicit VertexIDTrainingTreeMaker( const edm::ParameterSet & );
    ~VertexIDTrainingTreeMaker();

    static void fillDescriptions( edm::ConfigurationDescriptions &descriptions );


private:

    edm::Service<TFileService> fs_;

    virtual void beginJob() override;
    virtual void analyze( const edm::Event &, const edm::EventSetup & ) override;
    virtual void endJob() override;

    void initEventStructure();
    int mcTruthVertexIndex( const std::vector<edm::Ptr<reco::GenParticle> > &genParticles, const std::vector<edm::Ptr<reco::Vertex> > &vertices, const double &dzMatch = 0.1 );
    int sortedIndex( const unsigned int &trueVtxIndex, const unsigned int &sizemax, const Ptr<flashgg::DiPhotonCandidate> &diphoPtr );

    edm::EDGetTokenT<edm::View<flashgg::DiPhotonCandidate> > diphotonToken_;
    edm::EDGetTokenT<edm::View<reco::Vertex> >               vertexToken_;
    edm::EDGetTokenT<VertexCandidateMap> vertexCandidateMapTokenDz_;

    edm::EDGetTokenT<reco::BeamSpot > beamSpotToken_;
    edm::EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
    edm::EDGetTokenT<GenEventInfoProduct>  genEventInfoToken_;


    TTree *signalTree;
    TTree *backgroundTree;

    BackgroundInfo bkgInfo;
    SignalInfo sigInfo;

};

// ******************************************************************************************


//
// constructors and destructor
//
VertexIDTrainingTreeMaker::VertexIDTrainingTreeMaker( const edm::ParameterSet &iConfig ):
    diphotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
    vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
    vertexCandidateMapTokenDz_( consumes<VertexCandidateMap>( iConfig.getParameter<InputTag>( "VertexCandidateMapTagDz" ) ) ),
    beamSpotToken_( consumes<reco::BeamSpot >( iConfig.getParameter<InputTag>( "BeamSpotTag" ) ) ),
    genParticleToken_( consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag> ( "GenParticleTag" ) ) ),
    genEventInfoToken_( consumes<GenEventInfoProduct>( iConfig.getParameter<InputTag> ( "GenEventInfo" ) ) )
{
}

VertexIDTrainingTreeMaker::~VertexIDTrainingTreeMaker()
{
}

void
VertexIDTrainingTreeMaker::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup )
{

    // ********************************************************************************

    // access edm objects
    Handle<View<flashgg::DiPhotonCandidate> > diphotons;
    iEvent.getByToken( diphotonToken_, diphotons );
    const std::vector<edm::Ptr<flashgg::DiPhotonCandidate> > diphotonsPtrs = diphotons->ptrs();

    Handle<View<reco::Vertex> > primaryVertices;
    iEvent.getByToken( vertexToken_, primaryVertices );

    Handle<VertexCandidateMap> vertexCandidateMapDz;
    iEvent.getByToken( vertexCandidateMapTokenDz_, vertexCandidateMapDz );

    Handle<reco::BeamSpot> recoBeamSpotHandle;
    iEvent.getByToken( beamSpotToken_, recoBeamSpotHandle );
    math::XYZPoint beamSpot;
    if( recoBeamSpotHandle.isValid() ) {
        beamSpot = recoBeamSpotHandle->position();
    } else {
        cout << " WARNING NO VALID BEAM SPOT: this should not happen!" << endl;
    }

    Handle<View<reco::GenParticle> > genParticles;
    iEvent.getByToken( genParticleToken_, genParticles );
    const std::vector<edm::Ptr<reco::GenParticle> > genParticlesPtrs = genParticles->ptrs();

    Handle<GenEventInfoProduct> genEventInfo;
    iEvent.getByToken( genEventInfoToken_, genEventInfo );


    // ********************************************************************************

    initEventStructure();

    // diphoton loop

    for( size_t idipho = 0; idipho < diphotonsPtrs.size(); idipho++ ) {

        sigInfo.nvertex = primaryVertices->size();
        sigInfo.ndipho = diphotonsPtrs.size();
        bkgInfo.nvertex = primaryVertices->size();
        bkgInfo.ndipho = diphotonsPtrs.size();

        sigInfo.dipho_index = idipho;
        bkgInfo.dipho_index = idipho;
        sigInfo.genweight = genEventInfo->weight();
        bkgInfo.genweight = genEventInfo->weight();

        // get true vertex index:

        int trueVtxIndexI = mcTruthVertexIndex( genParticles->ptrs(), primaryVertices->ptrs() );
        if( trueVtxIndexI < 0 ) { continue; }

        Ptr<flashgg::DiPhotonCandidate> diphoPtr = diphotonsPtrs[idipho];
        sigInfo.dipho_mass = diphoPtr->mass();
        bkgInfo.dipho_mass = diphoPtr->mass();

        unsigned int trueVtxIndex = trueVtxIndexI;
        int trueVtxSortedIndexI = sortedIndex( trueVtxIndex, primaryVertices->size(), diphoPtr );
        if( trueVtxSortedIndexI < 0 ) { continue; }

        unsigned int trueVtxSortedIndex = trueVtxSortedIndexI;

        // Fill Signal Info
        if( trueVtxSortedIndex < diphoPtr->nVtxInfoSize() ) {
            sigInfo.LogSumPt2 = diphoPtr->logSumPt2( trueVtxSortedIndex );
            sigInfo.PtBal  =  diphoPtr->ptBal( trueVtxSortedIndex );
            sigInfo.PtAsym  =  diphoPtr->ptAsym( trueVtxSortedIndex );
            sigInfo.NConv  =  diphoPtr->nConv( trueVtxSortedIndex );
            sigInfo.PullConv  =  diphoPtr->pullConv( trueVtxSortedIndex );
            signalTree->Fill();
        }

        vector<int>	pvVecNoTrue;
        for( unsigned int i = 0 ; i < primaryVertices->size() ; i++ ) {
            if( i != trueVtxIndex ) { pvVecNoTrue.push_back( i ); }
        }

        int irand = -999;
        if( pvVecNoTrue.size() > 1 ) { irand = rand() % pvVecNoTrue.size(); }

        int randVtxIndex = -999;
        if( irand != -999 ) { randVtxIndex = pvVecNoTrue[irand]; }

        int randVtxSortedIndexI = sortedIndex( randVtxIndex, primaryVertices->size(), diphoPtr );
        if( randVtxSortedIndexI < 0 ) { continue; }
        unsigned int randVtxSortedIndex = randVtxSortedIndexI;

        // Fill Background Info

        if( randVtxSortedIndex < diphoPtr->nVtxInfoSize() ) {
            bkgInfo.LogSumPt2 = diphoPtr->logSumPt2( randVtxSortedIndex );
            bkgInfo.PtBal  =  diphoPtr->ptBal( randVtxSortedIndex );
            bkgInfo.PtAsym  =  diphoPtr->ptAsym( randVtxSortedIndex );
            bkgInfo.NConv  =  diphoPtr->nConv( randVtxSortedIndex );
            bkgInfo.PullConv  =  diphoPtr->pullConv( randVtxSortedIndex );
            backgroundTree->Fill();
        }

    }  // end diphoton candidate loop

}


void
VertexIDTrainingTreeMaker::beginJob()
{
    signalTree = fs_->make<TTree>( "signalTree", "per-diphoton tree" );
    signalTree->Branch( "nvertex"          , &sigInfo.nvertex         , "nvertex/I"        );
    signalTree->Branch( "ndipho"           , &sigInfo.ndipho          , "ndipho/I"         );
    signalTree->Branch( "genweight"        , &sigInfo.genweight       , "genweight/F"      );
    signalTree->Branch( "dipho_index"      , &sigInfo.dipho_index     , "dipho_index/I"    );
    signalTree->Branch( "dipho_mass"       , &sigInfo.dipho_mass      , "dipho_mass/F"     );
    signalTree->Branch( "logsumpt2"        , &sigInfo.LogSumPt2       , "logsumpt2/F"      );
    signalTree->Branch( "ptbal"            , &sigInfo.PtBal           , "ptbal/F"          );
    signalTree->Branch( "ptasym"           , &sigInfo.PtAsym          , "ptasym/F"         );
    signalTree->Branch( "nConv"            , &sigInfo.NConv           , "nConv/F"          );
    signalTree->Branch( "limPullToConv"    , &sigInfo.PullConv        , "limPullToConv/F"  );

    backgroundTree = fs_->make<TTree>( "backgroundTree", "per-diphoton tree" );
    backgroundTree->Branch( "nvertex"          , &bkgInfo.nvertex         , "nvertex/I"         );
    backgroundTree->Branch( "ndipho"           , &bkgInfo.ndipho          , "ndipho/I"          );
    backgroundTree->Branch( "genweight"        , &bkgInfo.genweight       , "genweight/F"       );
    backgroundTree->Branch( "dipho_index"      , &bkgInfo.dipho_index     , "dipho_index/I"     );
    backgroundTree->Branch( "dipho_mass"       , &bkgInfo.dipho_mass      , "dipho_mass/F"      );
    backgroundTree->Branch( "logsumpt2"        , &bkgInfo.LogSumPt2       , "logsumpt2/F"       );
    backgroundTree->Branch( "ptbal"            , &bkgInfo.PtBal           , "ptbal/F"           );
    backgroundTree->Branch( "ptasym"           , &bkgInfo.PtAsym          , "ptasym/F"          );
    backgroundTree->Branch( "nConv"            , &bkgInfo.NConv           , "nConv/F"           );
    backgroundTree->Branch( "limPullToConv"    , &bkgInfo.PullConv        , "limPullToConv/F"   );
}

void
VertexIDTrainingTreeMaker::endJob()
{
}

void
VertexIDTrainingTreeMaker::initEventStructure()
{

    sigInfo.nvertex         = -999;
    sigInfo.ndipho          = -999;
    sigInfo.dipho_index     = -999;
    sigInfo.dipho_mass      = -999.;

    sigInfo.LogSumPt2       = -999.;
    sigInfo.PtBal           = -999.;
    sigInfo.PtAsym          = -999.;
    sigInfo.NConv           = -999.;
    sigInfo.PullConv        = -999.;

    sigInfo.genweight       = -999.;

    bkgInfo.nvertex         = -999;
    bkgInfo.ndipho          = -999;
    bkgInfo.dipho_index     = -999;
    bkgInfo.dipho_mass      = -999.;

    bkgInfo.LogSumPt2       = -999.;
    bkgInfo.PtBal           = -999.;
    bkgInfo.PtAsym          = -999.;
    bkgInfo.NConv           = -999.;
    bkgInfo.PullConv        = -999.;

    bkgInfo.genweight       = -999.;
}

void
VertexIDTrainingTreeMaker::fillDescriptions( edm::ConfigurationDescriptions &descriptions )
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault( desc );
}

int 
VertexIDTrainingTreeMaker::sortedIndex( const unsigned int &trueVtxIndex, const unsigned int &sizemax, const Ptr<flashgg::DiPhotonCandidate> &diphoPtr )
{

    for( unsigned int j = 0; j < sizemax; j++ ) {
        int index = diphoPtr->mvaSortedIndex( j );
        if( index < 0 ) { continue; }
        if( ( unsigned int ) index == trueVtxIndex ) { return j; }
    }
    return -1;
}

int 
VertexIDTrainingTreeMaker::mcTruthVertexIndex( const std::vector<edm::Ptr<reco::GenParticle> > &genParticles ,
        const std::vector<edm::Ptr<reco::Vertex> > &vertices, const double &dzMatch )
{

    reco::Vertex::Point hardVertex( 0, 0, 0 );

    for( unsigned int genLoop = 0 ; genLoop < genParticles.size(); genLoop++ ) {

        if( fabs( genParticles[genLoop]->pdgId() ) < 10 || fabs( genParticles[genLoop]->pdgId() ) == 25 ) {
            hardVertex.SetCoordinates( genParticles[genLoop]->vx(), genParticles[genLoop]->vy(), genParticles[genLoop]->vz() );
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

DEFINE_FWK_MODULE( VertexIDTrainingTreeMaker );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

