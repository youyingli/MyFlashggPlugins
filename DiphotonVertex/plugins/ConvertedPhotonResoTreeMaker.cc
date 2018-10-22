
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

#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"

#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"

// **********************************************************************

// define the structures used to create tree branches and fill the trees

struct resolutionInfo {

    int ntracks;

    double dipho_mass;
    double pho_pt;
    int passEleVeto;
    double ztrue;
    double zvtx;
    double phoEtaSC;
    double gxy;
    double gz;
    double vtxZFromConv;
    double vtxZFromConvSC;

};

// **********************************************************************

using namespace std;
using namespace edm;
using namespace flashgg;


// **********************************************************************

class ConvertedPhotonResoTreeMaker : public edm::EDAnalyzer
{
public:
    explicit ConvertedPhotonResoTreeMaker( const edm::ParameterSet & );
    ~ConvertedPhotonResoTreeMaker();

    static void fillDescriptions( edm::ConfigurationDescriptions &descriptions );


private:

    edm::Service<TFileService> fs_;

    virtual void beginJob() override;
    virtual void analyze( const edm::Event &, const edm::EventSetup & ) override;
    virtual void endJob() override;

    void initEventStructure();
    std::vector<double> ZmcTruthVertexAndRecoVectex( const std::vector<edm::Ptr<reco::GenParticle> > &genParticles ,
            const std::vector<edm::Ptr<reco::Vertex> > &vertices, double dzMatch = 0.1);
    double vtxZFromConvOnly( const flashgg::Photon *pho, const edm::Ptr<reco:: Conversion> &conversion, const math::XYZPoint &beamSpot ) const;
    double vtxZFromConvSuperCluster( const flashgg::Photon *pho, const edm::Ptr<reco:: Conversion> &conversion, const math::XYZPoint &beamSpot ) const;
    std::vector<int> IndexMatchedConversion( const flashgg::Photon *g, 
            const std::vector<edm::Ptr<reco::Conversion> > &conversionsVector,  const std::vector<edm::Ptr<reco::Conversion> > &conversionsVectorSingleLeg,
            bool useSingleLeg ) const;



    edm::EDGetTokenT<edm::View<flashgg::DiPhotonCandidate> > diphotonToken_;
    edm::EDGetTokenT<edm::View<reco::Vertex> >          vertexToken_;
    edm::EDGetTokenT<View<reco::Conversion> >           conversionToken_;
    edm::EDGetTokenT<View<reco::Conversion> >           singlelegconversionToken_;
    edm::EDGetTokenT<reco::BeamSpot >                   beamSpotToken_;
    edm::EDGetTokenT<View<reco::GenParticle> >          genParticleToken_;

    TTree *resolutionTree;
    resolutionInfo resInfo;
};

// ******************************************************************************************


//
// constructors and destructor
//
ConvertedPhotonResoTreeMaker::ConvertedPhotonResoTreeMaker( const edm::ParameterSet &iConfig ):
    diphotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
    vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
    conversionToken_( consumes<View<reco::Conversion> >( iConfig.getParameter<InputTag>( "ConversionTag" ) ) ),
    singlelegconversionToken_( consumes<View<reco::Conversion> >( iConfig.getParameter<InputTag>( "SingleLegConversionTag" ) ) ),
    beamSpotToken_( consumes<reco::BeamSpot >( iConfig.getParameter<InputTag>( "BeamSpotTag" ) ) ),
    genParticleToken_( consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag> ( "GenParticleTag" ) ) )
{
}


ConvertedPhotonResoTreeMaker::~ConvertedPhotonResoTreeMaker()
{
}

void
ConvertedPhotonResoTreeMaker::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup )
{

    // ********************************************************************************

    // access edm objects
    Handle<View<flashgg::DiPhotonCandidate> > diphotons;
    iEvent.getByToken( diphotonToken_, diphotons );
    const std::vector<edm::Ptr<flashgg::DiPhotonCandidate> > diphotonsPtrs = diphotons->ptrs();

    Handle<View<reco::Vertex> > primaryVertices;
    iEvent.getByToken( vertexToken_, primaryVertices );

    Handle<View<reco::Conversion> > conversions;
    iEvent.getByToken(conversionToken_,conversions);

    Handle<View<reco::Conversion> > singlelegconversions;
    iEvent.getByToken(singlelegconversionToken_,singlelegconversions);
    
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


    // ********************************************************************************

    initEventStructure();

    vector<double> zposition = ZmcTruthVertexAndRecoVectex( genParticles->ptrs(), primaryVertices->ptrs() );

    if (diphotonsPtrs.size() > 0 && zposition.size() == 2) {
        resInfo.dipho_mass = diphotonsPtrs[0]->mass();
        resInfo.ztrue = zposition[0];
        resInfo.zvtx = zposition[1];
        vector<const flashgg::Photon*> diphoton;
        diphoton.push_back(diphotonsPtrs[0]->leadingPhoton());
        diphoton.push_back(diphotonsPtrs[0]->subLeadingPhoton());

        for (const auto& photon_it : diphoton) {
            const flashgg::Photon* photon = photon_it;
            resInfo.pho_pt = photon->pt();
            resInfo.passEleVeto = (int) photon->passElectronVeto();

            vector<int> indexConversion = IndexMatchedConversion(photon, conversions->ptrs(), singlelegconversions->ptrs(), true);
            if (indexConversion[1] == -1) continue;
            Ptr<reco:: Conversion> conversion = indexConversion[1] == 2? conversions->ptrAt(indexConversion[0]) : singlelegconversions->ptrAt(indexConversion[0]);
            double vtxZFromConv   = vtxZFromConvOnly(photon, conversion, beamSpot);
            double vtxZFromConvSC = vtxZFromConvSuperCluster(photon, conversion, beamSpot);

            resInfo.phoEtaSC = photon->superCluster()->eta();
            resInfo.gxy = sqrt( conversion->conversionVertex().x() * conversion->conversionVertex().x() + conversion->conversionVertex().y() *
                                conversion->conversionVertex().y() );
            resInfo.gz = conversion->conversionVertex().z();
            resInfo.vtxZFromConv = vtxZFromConv;
            resInfo.vtxZFromConvSC = vtxZFromConvSC;
            resInfo.ntracks = indexConversion[1];
            resolutionTree->Fill();
        }  // end photon candidate loop
    }
}


void
ConvertedPhotonResoTreeMaker::beginJob()
{
    resolutionTree = fs_->make<TTree>( "resolutionTree", "per-diphoton tree" );
    resolutionTree->Branch( "ntracks"          , &resInfo.ntracks,           "ntracks/I"          );
    resolutionTree->Branch( "dipho_mass"       , &resInfo.dipho_mass,        "dipho_mass/D"       );
    resolutionTree->Branch( "pho_pt"           , &resInfo.pho_pt,            "pho_pt/D"           );
    resolutionTree->Branch( "passEleVeto"      , &resInfo.passEleVeto,       "passEleVeto/I"      );
    resolutionTree->Branch( "ztrue"            , &resInfo.ztrue,             "ztrue/D"            );
    resolutionTree->Branch( "zvtx"             , &resInfo.zvtx,              "zvtx/D"             );
    resolutionTree->Branch( "phoEtaSC"         , &resInfo.phoEtaSC,          "phoEtaSC/D"         );
    resolutionTree->Branch( "gxy"              , &resInfo.gxy,               "gxy/D"              );
    resolutionTree->Branch( "gz"               , &resInfo.gz,                "gz/D"               );
    resolutionTree->Branch( "vtxZFromConv"     , &resInfo.vtxZFromConv,      "vtxZFromConv/D"     );
    resolutionTree->Branch( "vtxZFromConvSC"   , &resInfo.vtxZFromConvSC,    "vtxZFromConvSC/D"   );
}

void
ConvertedPhotonResoTreeMaker::endJob()
{
}

void
ConvertedPhotonResoTreeMaker::initEventStructure()
{
    resInfo.ntracks       = -999.;
    resInfo.dipho_mass    = -999.;
    resInfo.pho_pt        = -999.;
    resInfo.passEleVeto   = 0;
    resInfo.ztrue         = -999.;
    resInfo.zvtx          = -999.;
    resInfo.phoEtaSC      = -999.;
    resInfo.gxy           = -999.;
    resInfo.gz            = -999.;
    resInfo.vtxZFromConv  = -999.;
    resInfo.vtxZFromConvSC= -999.;
}

void
ConvertedPhotonResoTreeMaker::fillDescriptions( edm::ConfigurationDescriptions &descriptions )
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault( desc );
}

std::vector<double>
ConvertedPhotonResoTreeMaker::ZmcTruthVertexAndRecoVectex( const std::vector<edm::Ptr<reco::GenParticle> > &genParticles ,
        const std::vector<edm::Ptr<reco::Vertex> > &vertices, double dzMatch )
{

    reco::Vertex::Point hardVertex( 0, 0, 0 );

    for( unsigned int genLoop = 0 ; genLoop < genParticles.size(); genLoop++ ) {

        if( fabs( genParticles[genLoop]->pdgId() ) == 25 ) {
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

    std::vector<double> zposition;
    if( dzMin < dzMatch ) { 
        zposition.push_back(hardVertex.z());
        zposition.push_back(vertices[ivMatch]->z()); 
    }

    return zposition;
}

double 
//ConvertedPhotonResoTreeMaker::vtxZFromConvOnly( const edm::Ptr<flashgg::Photon> &pho, const edm::Ptr<reco:: Conversion> &conversion,
ConvertedPhotonResoTreeMaker::vtxZFromConvOnly( const flashgg::Photon *pho, const edm::Ptr<reco:: Conversion> &conversion,
        const math::XYZPoint &beamSpot ) const
{
    double dz = 0;
    if( conversion->nTracks() == 2 ) {
        double r = sqrt( conversion->refittedPairMomentum().perp2() );
        dz = ( conversion->conversionVertex().z() - beamSpot.z() )
             -
             ( ( conversion->conversionVertex().x() - beamSpot.x() ) * conversion->refittedPair4Momentum().x() + ( conversion->conversionVertex().y() - beamSpot.y() ) *
               conversion->refittedPair4Momentum().y() ) / r * conversion->refittedPair4Momentum().z() / r;
    }
    if( conversion->nTracks() == 1 ) {
        double r = sqrt( conversion->tracksPin()[0].x() * conversion->tracksPin()[0].x() + conversion->tracksPin()[0].y() * conversion->tracksPin()[0].y() );
        dz = ( conversion->conversionVertex().z() - beamSpot.z() )
             -
             ( ( conversion->conversionVertex().x() - beamSpot.x() ) * conversion->tracksPin()[0].x() + ( conversion->conversionVertex().y() - beamSpot.y() ) *
               conversion->tracksPin()[0].y() ) / r * conversion->tracksPin()[0].z() / r;
    }
    return dz + beamSpot.z();
}

double 
//ConvertedPhotonResoTreeMaker::vtxZFromConvSuperCluster( const edm::Ptr<flashgg::Photon> &pho, const edm::Ptr<reco:: Conversion> &conversion,
ConvertedPhotonResoTreeMaker::vtxZFromConvSuperCluster( const flashgg::Photon *pho, const edm::Ptr<reco:: Conversion> &conversion,
        const math::XYZPoint &beamSpot ) const
{
    // get the z from conversion plus SuperCluster
    double deltaX1 =  pho->caloPosition().x() - conversion->conversionVertex().x();
    double deltaY1 =  pho->caloPosition().y() - conversion->conversionVertex().y();
    double deltaZ1 =  pho->caloPosition().z() - conversion->conversionVertex().z();
    double R1 = sqrt( deltaX1 * deltaX1 + deltaY1 * deltaY1 );
    double tantheta = R1 / deltaZ1;
                                                                                                                                            
    double deltaX2 = conversion->conversionVertex().x() - beamSpot.x();
    double deltaY2 = conversion->conversionVertex().y() - beamSpot.y();
    double R2 = sqrt( deltaX2 * deltaX2 + deltaY2 * deltaY2 );
    double deltaZ2 = R2 / tantheta;
    double higgsZ =  pho->caloPosition().z() - deltaZ1 - deltaZ2;
    return higgsZ;
}

std::vector<int> 
//ConvertedPhotonResoTreeMaker::IndexMatchedConversion( const edm::Ptr<flashgg::Photon> &g,
ConvertedPhotonResoTreeMaker::IndexMatchedConversion( const flashgg::Photon *g,
        const std::vector<edm::Ptr<reco::Conversion> > &conversionsVector,  const std::vector<edm::Ptr<reco::Conversion> > &conversionsVectorSingleLeg, 
        bool useSingleLeg ) const
{
    double mindR = 999;
    int nConvLegs = 0;
    bool doOneLeg = true;
    bool pureGeomConvMatching = true;
                                                                                                                                                        
    std::vector<int> result;
                                                                                                                                                        
    if(!pureGeomConvMatching) assert( g->hasConversionTracks() );
    int selected_conversion_index = -1;
                                                                                                                                                        
    if( (g->hasConversionTracks() && !pureGeomConvMatching) || pureGeomConvMatching){
        
        for( unsigned int i = 0; i < conversionsVector.size(); i++ ) {
            edm::Ptr<reco::Conversion> conv = conversionsVector[i];
            if( conv->nTracks() == 2 ) {
                if( !conv->isConverted() ) { continue; }
                if( conv->refittedPair4Momentum().pt() < 10. ) { continue; }
                if( TMath::Prob( conv->conversionVertex().chi2(), conv->conversionVertex().ndof() ) < 1e-6 ) { continue; }
                
                TVector3 VtxtoSC;
                VtxtoSC.SetXYZ( g->superCluster()->position().x() - conv->conversionVertex().x(),
                                g->superCluster()->position().y() - conv->conversionVertex().y(),
                                g->superCluster()->position().z() - conv->conversionVertex().z() );
                TVector3 RefPairMo;
                RefPairMo.SetXYZ( conv->refittedPairMomentum().x(), conv->refittedPairMomentum().y(), conv->refittedPairMomentum().z() );
                double dR = 0;
                dR = VtxtoSC.DeltaR( RefPairMo );
                if( dR < mindR ) {
                    mindR = dR;
                    selected_conversion_index = i;
                }
            }
        }
        if( mindR < 0.1 ) {
            result.push_back( selected_conversion_index );
            nConvLegs = 2;
            result.push_back( nConvLegs );
            doOneLeg = false;
        }
        if( doOneLeg && useSingleLeg ) {
            mindR = 999;
            for( unsigned int j = 0; j < conversionsVectorSingleLeg.size(); j++ ) {
                edm::Ptr<reco::Conversion> conv = conversionsVectorSingleLeg[j];
                if( conv->nTracks() == 1 ) {
                    TVector3 VtxtoSC;
                    VtxtoSC.SetXYZ( g->superCluster()->position().x() - conv->conversionVertex().x(),
                                    g->superCluster()->position().y() - conv->conversionVertex().y(),
                                    g->superCluster()->position().z() - conv->conversionVertex().z() );
                    TVector3 RefPairMo;
                    float oneLegTrack_X = conv->tracksPin()[0].x();
                    float oneLegTrack_Y = conv->tracksPin()[0].y();
                    float oneLegTrack_Z = conv->tracksPin()[0].z();
                    
                    RefPairMo.SetXYZ( oneLegTrack_X, oneLegTrack_Y, oneLegTrack_Z );
                    double dR = 0;
                    dR = VtxtoSC.DeltaR( RefPairMo );
                    if( dR < mindR ) {
                        mindR = dR;
                        selected_conversion_index = j;
                    }                        
                }
            }
            if( mindR < 0.1 ) {
                result.push_back( selected_conversion_index );
                nConvLegs = 1;
                result.push_back( nConvLegs );
            }
        }
    }
    
    if( mindR < 0.1 )
        {return result;}
    else {
        result.push_back( -1 );
        result.push_back( -1 );
        return result;
    }
}

DEFINE_FWK_MODULE( ConvertedPhotonResoTreeMaker );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

