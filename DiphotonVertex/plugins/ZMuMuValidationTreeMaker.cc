
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

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
//for miniAOD:
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "flashgg/MicroAOD/interface/VertexSelectorBase.h"
#include "flashgg/DataFormats/interface/VertexCandidateMap.h"

#include "TTree.h"

// **********************************************************************

// define the structures used to create tree branches and fill the trees


// **********************************************************************

using namespace std;
using namespace edm;
using namespace flashgg;

struct GlobalInfo {
    int nPU               ;
    float BSsigmaz        ;
};

struct DimuonInfo {

    float BSz             ;
    float BSsigmaz        ;
    float mcWeight        ;
    float zRecoTrue       ;    
    float zRecoTrueNoMu   ;
    float zTrue           ;
    int   nPU             ;
    int   ntightmuons     ;
    float dimuon_mass     ; 
    float dimuon_pt       ; 
    float muon1_pt        ; 
    float muon1_eta       ; 
    float muon2_pt        ; 
    float muon2_eta       ; 
    float zChosenWithMu   ; 
    float zChosenNoMu     ; 
    float zRecoFromMu     ;
    float zRecoFromMuNoMu ;

    float mvaProbWithMu   ; 
    float logSumPt2WithMu ; 
    float ptBalWithMu     ; 
    float ptAsymWithMu    ; 
    int   nVtxWithMu      ; 
    float mva0WithMu      ; 
    float mva1WithMu      ; 
    float mva2WithMu      ; 
    float dZ1WithMu       ; 
    float dZ2WithMu       ;

    float mvaProbNoMu     ; 
    float logSumPt2NoMu   ; 
    float ptBalNoMu       ; 
    float ptAsymNoMu      ; 
    int   nVtxNoMu        ; 
    float mva0NoMu        ; 
    float mva1NoMu        ; 
    float mva2NoMu        ; 
    float dZ1NoMu         ; 
    float dZ2NoMu         ; 

};

struct SignalInfo {

    float vtxIdmva;
    float LogSumPt2;
    float PtBal;
    float PtAsym;
    float DZtrue;

};

struct BackgroundInfo {

    float vtxIdmva;
    float LogSumPt2;
    float PtBal;
    float PtAsym;
    float DZtrue;

};

// **********************************************************************

class ZMuMuValidationTreeMaker : public edm::EDAnalyzer
{
public:
    explicit ZMuMuValidationTreeMaker( const edm::ParameterSet & );
    ~ZMuMuValidationTreeMaker();
    
    static void fillDescriptions( edm::ConfigurationDescriptions &descriptions );
    
    
private:
    
    edm::Service<TFileService> fs_;
    
    virtual void beginJob() override;
    virtual void analyze( const edm::Event &, const edm::EventSetup & ) override;
    virtual void endJob() override;

    void initEventStructure();
   
    float getMCTruthVertexZ( const std::vector<edm::Ptr<reco::GenParticle>> &gens );
    int getRecoClosestToTrueVertexIndex(  const std::vector<edm::Ptr<reco::GenParticle>> &gens, const vector<edm::Ptr<reco::Vertex> > &vertices, const double &dzMatch = 0.2 );
    int getRecoWithMuonsVertexIndex( const vector<edm::Ptr<reco::Vertex> > &vertices, const Ptr<pat::Muon> &mu1, const Ptr<pat::Muon> &mu2, const double &dzMatch = 0.2 );
    int sortedIndex( const unsigned int &trueVtxIndex, const std::vector<int> &VtxSortedIndexVector);
    std::vector<int> isTightMuonWithoutVtxCut(const vector<edm::Ptr<pat::Muon>> &patMuons); 

    edm::EDGetTokenT<View<pat::Muon> > patMuonToken_;
    edm::EDGetTokenT<edm::View<reco::Vertex> >               vertexTokenWithMu_;
    edm::EDGetTokenT<edm::View<reco::Vertex> >               vertexTokenNoMu_;
    edm::EDGetTokenT< VertexCandidateMap > vertexCandidateMapToken_;
    edm::EDGetTokenT< VertexCandidateMap > vertexCandidateMapNoMuToken_;
    unique_ptr<VertexSelectorBase> vertexSelector_;
    edm::EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
    edm::EDGetToken genEventInfoToken_;
    edm::EDGetTokenT<reco::BeamSpot > beamSpotToken_;
    edm::EDGetTokenT<edm::View<PileupSummaryInfo> >  PileUpToken_;

    TTree *GlobalTree;
    TTree *DiMuonTree;
    TTree *SignalTree;
    TTree *BackgroundTree;

    GlobalInfo globalInfo;
    DimuonInfo dimuonInfo;
    SignalInfo sigInfo;
    BackgroundInfo bkgInfo;
};

// ******************************************************************************************


//
// constructors and destructor
//
ZMuMuValidationTreeMaker::ZMuMuValidationTreeMaker( const edm::ParameterSet &iConfig ):
    patMuonToken_( consumes<View<pat::Muon> >( iConfig.getParameter<InputTag>( "MuonTag" ) ) ),
    vertexTokenWithMu_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTagWithMu" ) ) ),
    vertexTokenNoMu_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTagNoMu" ) ) ),
    vertexCandidateMapToken_( consumes<VertexCandidateMap>( iConfig.getParameter<InputTag>( "VertexCandidateMapTag" ) ) ),
    vertexCandidateMapNoMuToken_( consumes<VertexCandidateMap>( iConfig.getParameter<InputTag>( "VertexCandidateMapTagNoMu" ) ) ),
    genParticleToken_( consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag> ( "GenParticleTag" ) ) ),
    genEventInfoToken_( consumes<GenEventInfoProduct>( iConfig.getParameter <edm::InputTag> ( "genEventInfoProduct" ) ) ),
    beamSpotToken_( consumes<reco::BeamSpot >( iConfig.getParameter<InputTag>( "BeamSpotTag" ) ) ),
    PileUpToken_( consumes<View<PileupSummaryInfo> >( iConfig.getParameter<InputTag> ( "PileUpTag" ) ) )
{
        const std::string &VertexSelectorName = iConfig.getParameter<std::string>( "VertexSelectorName" );
        vertexSelector_.reset( FlashggVertexSelectorFactory::get()->create( VertexSelectorName, iConfig ) );
}


ZMuMuValidationTreeMaker::~ZMuMuValidationTreeMaker()
{


}



void
ZMuMuValidationTreeMaker::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup )
{
    Handle<View<reco::GenParticle> > genParticles;
    iEvent.getByToken( genParticleToken_, genParticles );

    Handle<GenEventInfoProduct> genEventInfo;
    iEvent.getByToken( genEventInfoToken_, genEventInfo );

    Handle<View<pat::Muon> >  patMuons;
    iEvent.getByToken( patMuonToken_, patMuons );

    Handle<View<reco::Vertex> > primaryVerticesNoMu;
    iEvent.getByToken( vertexTokenNoMu_, primaryVerticesNoMu );

    Handle<View<reco::Vertex> > primaryVerticesWithMu;
    iEvent.getByToken( vertexTokenWithMu_, primaryVerticesWithMu );

    Handle<VertexCandidateMap> vertexCandidateMap;
    iEvent.getByToken( vertexCandidateMapToken_, vertexCandidateMap );

    Handle<VertexCandidateMap> vertexCandidateMapNoMu;
    iEvent.getByToken( vertexCandidateMapNoMuToken_, vertexCandidateMapNoMu );

    Handle<reco::BeamSpot> recoBeamSpotHandle;
    iEvent.getByToken( beamSpotToken_, recoBeamSpotHandle );

    initEventStructure();

    math::XYZPoint beamSpot;
    double BSsigmaz = 0;
    if( recoBeamSpotHandle.isValid() ) {
        beamSpot = recoBeamSpotHandle->position();
        BSsigmaz = recoBeamSpotHandle->sigmaZ();
        globalInfo.BSsigmaz = BSsigmaz;
    } else {
        cout << " WARNING NO VALID BEAM SPOT: this should not happen!" << endl;
    }

    Handle<View< PileupSummaryInfo> > PileupInfos;
    iEvent.getByToken( PileUpToken_, PileupInfos ); 

    if( !iEvent.isRealData() ) {

        dimuonInfo.mcWeight = genEventInfo->weight();

        int trueRecoVtxIdx = getRecoClosestToTrueVertexIndex( genParticles->ptrs(), primaryVerticesWithMu->ptrs() );
        int trueRecoVtxIdxNoMu = getRecoClosestToTrueVertexIndex( genParticles->ptrs(), primaryVerticesNoMu->ptrs() );

        if(trueRecoVtxIdx != -1)     dimuonInfo.zRecoTrue     = primaryVerticesWithMu->ptrAt( trueRecoVtxIdx )->z();    
        if(trueRecoVtxIdxNoMu != -1) dimuonInfo.zRecoTrueNoMu = primaryVerticesNoMu->ptrAt( trueRecoVtxIdxNoMu )->z();
        dimuonInfo.zTrue = getMCTruthVertexZ( genParticles->ptrs() );
        
        // pileup info
        for( unsigned int PVI = 0; PVI < PileupInfos->size(); ++PVI ) {
            Int_t pu_bunchcrossing = PileupInfos->ptrAt( PVI )->getBunchCrossing();
            if( pu_bunchcrossing == 0 ) {
                globalInfo.nPU = PileupInfos->ptrAt( PVI )->getTrueNumInteractions();
                dimuonInfo.nPU = PileupInfos->ptrAt( PVI )->getTrueNumInteractions();
                break;
            }
        }
    }
    GlobalTree->Fill();
    
    
    vector<int> muon_index = isTightMuonWithoutVtxCut(patMuons->ptrs());
  
    if( muon_index.size() >= 2 ){
        
        float Zmass = 91.19;
        Ptr<pat::Muon> pat_muon1 = patMuons->ptrAt(muon_index[0]);
        Ptr<pat::Muon> pat_muon2 = patMuons->ptrAt(muon_index[1]);
   
        math::XYZTLorentzVector dimuonP4 = pat_muon1->p4() + pat_muon2->p4();
        double dimuon_mass = dimuonP4.M();
        double dimuon_pt = dimuonP4.pt();

        if ( fabs(dimuon_mass - Zmass) < 30. ) {
            dimuonInfo.ntightmuons = muon_index.size();
            dimuonInfo.dimuon_mass = dimuon_mass;
            dimuonInfo.dimuon_pt = dimuon_pt;

            dimuonInfo.BSz = beamSpot.z();
            dimuonInfo.BSsigmaz = BSsigmaz;

            dimuonInfo.muon1_pt  = pat_muon1->pt();
            dimuonInfo.muon1_eta = pat_muon1->eta();
            dimuonInfo.muon2_pt  = pat_muon2->pt();
            dimuonInfo.muon2_eta = pat_muon2->eta();

            //Selected Vertex
            vector<float> infoWithMu;
            vector<float> infoNoMu;

            Ptr<reco::Vertex> selectedVtxWithMu = vertexSelector_->select( pat_muon1, pat_muon2, primaryVerticesWithMu->ptrs(), *vertexCandidateMap, beamSpot ); 
            vertexSelector_->getInfoFromLastSelection(infoWithMu);

            Ptr<reco::Vertex> selectedVtxNoMu = vertexSelector_->select( pat_muon1, pat_muon2, primaryVerticesNoMu->ptrs(), *vertexCandidateMapNoMu, beamSpot ); 
            vertexSelector_->getInfoFromLastSelection(infoNoMu);
            
            dimuonInfo.zChosenWithMu = selectedVtxWithMu->z();
            dimuonInfo.zChosenNoMu = selectedVtxNoMu->z();

            dimuonInfo.mvaProbWithMu   = infoWithMu[0]; 
            dimuonInfo.logSumPt2WithMu = infoWithMu[1];
            dimuonInfo.ptBalWithMu     = infoWithMu[2];
            dimuonInfo.ptAsymWithMu    = infoWithMu[3];
            dimuonInfo.nVtxWithMu      = infoWithMu[4];
            dimuonInfo.mva0WithMu      = infoWithMu[5];
            dimuonInfo.mva1WithMu      = infoWithMu[6];
            dimuonInfo.mva2WithMu      = infoWithMu[7];
            dimuonInfo.dZ1WithMu       = infoWithMu[8];
            dimuonInfo.dZ2WithMu       = infoWithMu[9];

            dimuonInfo.mvaProbNoMu     = infoNoMu[0]; 
            dimuonInfo.logSumPt2NoMu   = infoNoMu[1];
            dimuonInfo.ptBalNoMu       = infoNoMu[2];
            dimuonInfo.ptAsymNoMu      = infoNoMu[3];
            dimuonInfo.nVtxNoMu        = infoNoMu[4];
            dimuonInfo.mva0NoMu        = infoNoMu[5];
            dimuonInfo.mva1NoMu        = infoNoMu[6];
            dimuonInfo.mva2NoMu        = infoNoMu[7];
            dimuonInfo.dZ1NoMu         = infoNoMu[8];
            dimuonInfo.dZ2NoMu         = infoNoMu[9];


            //Reconstruction real vertex
            int iFromMuonsWithMu = getRecoWithMuonsVertexIndex(primaryVerticesWithMu->ptrs(), pat_muon1, pat_muon2);
            int iFromMuonsNoMu = getRecoWithMuonsVertexIndex(primaryVerticesNoMu->ptrs(), pat_muon1, pat_muon2);       
      
            if( iFromMuonsWithMu != -1 ) {

                dimuonInfo.zRecoFromMu = primaryVerticesWithMu->ptrAt(iFromMuonsWithMu)->z();
                if( iFromMuonsNoMu != -1 ) dimuonInfo.zRecoFromMuNoMu = primaryVerticesNoMu->ptrAt(iFromMuonsNoMu)->z();
                DiMuonTree->Fill();

                //Right Vertex
                vector<int> VtxSortedIndexVector;
                for( int i = 0; i < (int)primaryVerticesNoMu->size(); i++) VtxSortedIndexVector.push_back( vertexSelector_->getSortedIndexFromLastSelection( i ) );
                unsigned int trueVtxIndex = iFromMuonsNoMu;
                int trueVtxSortedIndex = sortedIndex(trueVtxIndex, VtxSortedIndexVector);

                if (trueVtxSortedIndex != -1) {    
                    vector<float> info;
                    vertexSelector_->getInfoFromLastSelectionForVtxIdx( info , (unsigned int) trueVtxSortedIndex );
                    sigInfo.vtxIdmva  = info[0];
                    sigInfo.LogSumPt2 = info[1];
                    sigInfo.PtBal     = info[2];
                    sigInfo.PtAsym    = info[3];
                    sigInfo.DZtrue    = primaryVerticesNoMu->ptrAt(trueVtxIndex)->position().z() - primaryVerticesWithMu->ptrAt(iFromMuonsWithMu)->position().z();

                    SignalTree->Fill();
                }

                //Wrong Vertex
                vector<int>	notrueVtxIndexVector;
                for( unsigned int i = 0 ; i < primaryVerticesNoMu->size() ; i++ ) {
                    if( i != trueVtxIndex ) { notrueVtxIndexVector.push_back( i ); }
                }

                int irand = -999;
                if( notrueVtxIndexVector.size() > 1 ) { irand = rand() % notrueVtxIndexVector.size(); }

                int randVtxIndex = -999;
                if( irand != -999 ) { randVtxIndex = notrueVtxIndexVector[irand]; }

                int randVtxSortedIndex = sortedIndex(randVtxIndex, VtxSortedIndexVector);

                if( randVtxSortedIndex != -1 ) {
                    vector<float> info;
                    vertexSelector_->getInfoFromLastSelectionForVtxIdx( info , (unsigned int) randVtxSortedIndex );
                    bkgInfo.vtxIdmva  = info[0];
                    bkgInfo.LogSumPt2 = info[1];
                    bkgInfo.PtBal     = info[2];
                    bkgInfo.PtAsym    = info[3];
                    bkgInfo.DZtrue    = primaryVerticesNoMu->ptrAt(randVtxIndex)->position().z() - primaryVerticesWithMu->ptrAt(iFromMuonsWithMu)->position().z();

                    BackgroundTree->Fill();
                }
            }

        }//mass cut
    }//tight muon

}


void
ZMuMuValidationTreeMaker::beginJob()
{

    GlobalTree = fs_->make<TTree>( "GlobalTree", "GlobalTree" );   
    GlobalTree->Branch( "nPU"              , &globalInfo.nPU              , "nPU/I"             );  
    GlobalTree->Branch( "BSsigmaz"         , &globalInfo.BSsigmaz         , "BSsigmaz/F"        );  

    DiMuonTree = fs_->make<TTree>( "DiMuonTree", "DiMuonTree" );   
    DiMuonTree->Branch( "BSz"              , &dimuonInfo.BSz              , "BSz/F"             );  
    DiMuonTree->Branch( "BSsigmaz"         , &dimuonInfo.BSsigmaz         , "BSsigmaz/F"        );  
    DiMuonTree->Branch( "mcWeight"         , &dimuonInfo.mcWeight         , "mcWeight/F"        );  
    DiMuonTree->Branch( "zRecoTrue"        , &dimuonInfo.zRecoTrue        , "zRecoTrue/F"       );  
    DiMuonTree->Branch( "zRecoTrueNoMu"    , &dimuonInfo.zRecoTrueNoMu    , "zRecoTrueNoMu/F"   );  
    DiMuonTree->Branch( "zTrue"            , &dimuonInfo.zTrue            , "zTrue/F"           );  
    DiMuonTree->Branch( "nPU"              , &dimuonInfo.nPU              , "nPU/I"             );  
    DiMuonTree->Branch( "ntightmuons"      , &dimuonInfo.ntightmuons      , "ntightmuons/I"     );  
    DiMuonTree->Branch( "dimuon_mass"      , &dimuonInfo.dimuon_mass      , "dimuon_mass/F"     );  
    DiMuonTree->Branch( "dimuon_pt"        , &dimuonInfo.dimuon_pt        , "dimuon_pt/F"       );  
    DiMuonTree->Branch( "muon1_pt"         , &dimuonInfo.muon1_pt         , "muon1_pt/F"        );  
    DiMuonTree->Branch( "muon1_eta"        , &dimuonInfo.muon1_eta        , "muon1_eta/F"       );  
    DiMuonTree->Branch( "muon2_pt"         , &dimuonInfo.muon2_pt         , "muon2_pt/F"        );  
    DiMuonTree->Branch( "muon2_eta"        , &dimuonInfo.muon2_eta        , "muon2_eta/F"       );  
    DiMuonTree->Branch( "zChosenWithMu"    , &dimuonInfo.zChosenWithMu    , "zChosenWithMu/F"   );  
    DiMuonTree->Branch( "zChosenNoMu"      , &dimuonInfo.zChosenNoMu      , "zChosenNoMu/F"     );  
    DiMuonTree->Branch( "zRecoFromMu"      , &dimuonInfo.zRecoFromMu      , "zRecoFromMu/F"     );  
    DiMuonTree->Branch( "zRecoFromMuNoMu"  , &dimuonInfo.zRecoFromMuNoMu  , "zRecoFromMuNoMu/F" );  
    DiMuonTree->Branch( "mvaProbWithMu"    , &dimuonInfo.mvaProbWithMu    , "mvaProbWithMu/F"   );  
    DiMuonTree->Branch( "logSumPt2WithMu"  , &dimuonInfo.logSumPt2WithMu  , "logSumPt2WithMu/F" );  
    DiMuonTree->Branch( "ptBalWithMu"      , &dimuonInfo.ptBalWithMu      , "ptBalWithMu/F"     );  
    DiMuonTree->Branch( "ptAsymWithMu"     , &dimuonInfo.ptAsymWithMu     , "ptAsymWithMu/F"    );  
    DiMuonTree->Branch( "nVtxWithMu"       , &dimuonInfo.nVtxWithMu       , "nVtxWithMu/I"      );  
    DiMuonTree->Branch( "mva0WithMu"       , &dimuonInfo.mva0WithMu       , "mva0WithMu/F"      );  
    DiMuonTree->Branch( "mva1WithMu"       , &dimuonInfo.mva1WithMu       , "mva1WithMu/F"      );  
    DiMuonTree->Branch( "mva2WithMu"       , &dimuonInfo.mva2WithMu       , "mva2WithMu/F"      );  
    DiMuonTree->Branch( "dZ1WithMu"        , &dimuonInfo.dZ1WithMu        , "dZ1WithMu/F"       );  
    DiMuonTree->Branch( "dZ2WithMu"        , &dimuonInfo.dZ2WithMu        , "dZ2WithMu/F"       );  
    DiMuonTree->Branch( "mvaProbNoMu"      , &dimuonInfo.mvaProbNoMu      , "mvaProbNoMu/F"     );  
    DiMuonTree->Branch( "logSumPt2NoMu"    , &dimuonInfo.logSumPt2NoMu    , "logSumPt2NoMu/F"   );  
    DiMuonTree->Branch( "ptBalNoMu"        , &dimuonInfo.ptBalNoMu        , "ptBalNoMu/F"       );  
    DiMuonTree->Branch( "ptAsymNoMu"       , &dimuonInfo.ptAsymNoMu       , "ptAsymNoMu/F"      );  
    DiMuonTree->Branch( "nVtxNoMu"         , &dimuonInfo.nVtxNoMu         , "nVtxNoMu/I"        );  
    DiMuonTree->Branch( "mva0NoMu"         , &dimuonInfo.mva0NoMu         , "mva0NoMu/F"        );  
    DiMuonTree->Branch( "mva1NoMu"         , &dimuonInfo.mva1NoMu         , "mva1NoMu/F"        );  
    DiMuonTree->Branch( "mva2NoMu"         , &dimuonInfo.mva2NoMu         , "mva2NoMu/F"        );  
    DiMuonTree->Branch( "dZ1NoMu"          , &dimuonInfo.dZ1NoMu          , "dZ1NoMu/F"         );  
    DiMuonTree->Branch( "dZ2NoMu"          , &dimuonInfo.dZ2NoMu          , "dZ2NoMu/F"         );  

    SignalTree = fs_->make<TTree>( "SignalTree", "SignalTree" );
    SignalTree->Branch( "BSz"        , &dimuonInfo.BSz       , "BSz/F"      );
    SignalTree->Branch( "BSsigmaz"   , &dimuonInfo.BSsigmaz  , "BSsigmaz/F" );
    SignalTree->Branch( "mcWeight"   , &dimuonInfo.mcWeight  , "mcWeight/F" );
    SignalTree->Branch( "nPU"        , &dimuonInfo.nPU       , "nPU/I"      );
    SignalTree->Branch( "nvertex"    , &dimuonInfo.nVtxNoMu  , "nvertex/I"  );
    SignalTree->Branch( "vtxIdmva"   , &sigInfo.vtxIdmva     , "vtxIdmva/F" );
    SignalTree->Branch( "LogSumPt2"  , &sigInfo.LogSumPt2    , "LogSumPt2/F");
    SignalTree->Branch( "PtBal"      , &sigInfo.PtBal        , "PtBal/F"    );
    SignalTree->Branch( "PtAsym"     , &sigInfo.PtAsym       , "PtAsym/F"   );
    SignalTree->Branch( "DZtrue"     , &sigInfo.DZtrue       , "DZtrue/F"   );

    BackgroundTree = fs_->make<TTree>( "backgroundTree", "BackgroundTree" );
    BackgroundTree->Branch( "BSz"        , &dimuonInfo.BSz       , "BSz/F"      );
    BackgroundTree->Branch( "BSsigmaz"   , &dimuonInfo.BSsigmaz  , "BSsigmaz/F" );
    BackgroundTree->Branch( "mcWeight"   , &dimuonInfo.mcWeight  , "mcWeight/F" );
    BackgroundTree->Branch( "nPU"        , &dimuonInfo.nPU       , "nPU/I"      );
    BackgroundTree->Branch( "nvertex"    , &dimuonInfo.nVtxNoMu  , "nvertex/I"  );
    BackgroundTree->Branch( "vtxIdmva"   , &bkgInfo.vtxIdmva     , "vtxIdmva/F" );
    BackgroundTree->Branch( "LogSumPt2"  , &bkgInfo.LogSumPt2    , "LogSumPt2/F");
    BackgroundTree->Branch( "PtBal"      , &bkgInfo.PtBal        , "PtBal/F"    );
    BackgroundTree->Branch( "PtAsym"     , &bkgInfo.PtAsym       , "PtAsym/F"   );
    BackgroundTree->Branch( "DZtrue"     , &bkgInfo.DZtrue       , "DZtrue/F"   );

}

void
ZMuMuValidationTreeMaker::endJob()
{
}

void
ZMuMuValidationTreeMaker::initEventStructure()
{
    globalInfo.nPU             = -999;
    globalInfo.BSsigmaz        = -999.; 

    dimuonInfo.BSz             = -999.; 
    dimuonInfo.BSsigmaz        = -999.; 
    dimuonInfo.mcWeight        = -999.; 
    dimuonInfo.zRecoTrue       = -999.; 
    dimuonInfo.zRecoTrueNoMu   = -999.; 
    dimuonInfo.zTrue           = -999.; 
    dimuonInfo.nPU             = -999 ; 
    dimuonInfo.ntightmuons     = -999 ; 
    dimuonInfo.dimuon_mass     = -999.; 
    dimuonInfo.dimuon_pt       = -999.; 
    dimuonInfo.muon1_pt        = -999.; 
    dimuonInfo.muon1_eta       = -999.; 
    dimuonInfo.muon2_pt        = -999.; 
    dimuonInfo.muon2_eta       = -999.; 
    dimuonInfo.zChosenWithMu   = -999.; 
    dimuonInfo.zChosenNoMu     = -999.; 
    dimuonInfo.zRecoFromMu     = -999.; 
    dimuonInfo.zRecoFromMuNoMu = -999.; 
    dimuonInfo.mvaProbWithMu   = -999.; 
    dimuonInfo.logSumPt2WithMu = -999.; 
    dimuonInfo.ptBalWithMu     = -999.; 
    dimuonInfo.ptAsymWithMu    = -999.; 
    dimuonInfo.nVtxWithMu      = -999 ; 
    dimuonInfo.mva0WithMu      = -999.; 
    dimuonInfo.mva1WithMu      = -999.; 
    dimuonInfo.mva2WithMu      = -999.; 
    dimuonInfo.dZ1WithMu       = -999.; 
    dimuonInfo.dZ2WithMu       = -999.; 
    dimuonInfo.mvaProbNoMu     = -999.; 
    dimuonInfo.logSumPt2NoMu   = -999.; 
    dimuonInfo.ptBalNoMu       = -999.; 
    dimuonInfo.ptAsymNoMu      = -999.; 
    dimuonInfo.nVtxNoMu        = -999 ; 
    dimuonInfo.mva0NoMu        = -999.; 
    dimuonInfo.mva1NoMu        = -999.; 
    dimuonInfo.mva2NoMu        = -999.; 
    dimuonInfo.dZ1NoMu         = -999.; 
    dimuonInfo.dZ2NoMu         = -999.; 

    sigInfo.vtxIdmva   = -999.; 
    sigInfo.LogSumPt2  = -999.; 
    sigInfo.PtBal      = -999.; 
    sigInfo.PtAsym     = -999.; 
    sigInfo.DZtrue     = -999.; 

    bkgInfo.vtxIdmva   = -999.; 
    bkgInfo.LogSumPt2  = -999.; 
    bkgInfo.PtBal      = -999.; 
    bkgInfo.PtAsym     = -999.; 
    bkgInfo.DZtrue     = -999.; 

}

int 
ZMuMuValidationTreeMaker::getRecoWithMuonsVertexIndex( const vector<edm::Ptr<reco::Vertex> > &vertices, const Ptr<pat::Muon> &mu1, const Ptr<pat::Muon> &mu2, const double &dzMatch)
{

    double IP1 = mu1->vz();
    double IP2 = mu2->vz(); 

    double average = 0.5 * (IP1 + IP2);
        
    int  ivMatch = 0;
    int  ivMatch1 = 0;
    int  ivMatch2 = 0;
    double dzMin = 999;
    double dzMin1 = 999;
    double dzMin2 = 999;
    
    for( unsigned int iv = 0; iv < vertices.size(); iv++ ) {
        double dz = fabs( vertices[iv]->z() - average );
        double dz1 = fabs( vertices[iv]->z() - IP1 );
        double dz2 = fabs( vertices[iv]->z() - IP2 );
        if( dz < dzMin ) {
            ivMatch = iv;
            dzMin   = dz;
        } 
        if( dz1 < dzMin1 ) {
            ivMatch1 = iv;
            dzMin1   = dz1;
        }
        if( dz2 < dzMin2 ) {
            ivMatch2 = iv;
            dzMin2   = dz2;
        }
    }
     
    if(ivMatch == ivMatch1 && ivMatch == ivMatch2 && dzMin < dzMatch) return ivMatch;
    return -1;
}


float 
ZMuMuValidationTreeMaker::getMCTruthVertexZ( const std::vector<edm::Ptr<reco::GenParticle>> &gens )
{
    float zTrue = 999.;
    for( unsigned int genLoop = 0 ; genLoop < gens.size(); genLoop++ ) {
        if( fabs( gens[genLoop]->pdgId() ) == 23 ) {
            zTrue = gens[genLoop]->vz();
            break;
        }
    }
    return zTrue;
}

int ZMuMuValidationTreeMaker::getRecoClosestToTrueVertexIndex( const std::vector<edm::Ptr<reco::GenParticle>> &gens , const vector<edm::Ptr<reco::Vertex> > &vertices, const double &dzMatch )
{

    reco::Vertex::Point hardVertex( 0, 0, 0 );

    for( unsigned int genLoop = 0 ; genLoop < gens.size(); genLoop++ ) {

        if( fabs( gens[genLoop]->pdgId() ) < 10 || fabs( gens[genLoop]->pdgId() ) == 23 ) {
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

int
ZMuMuValidationTreeMaker::sortedIndex( const unsigned int &trueVtxIndex, const vector<int> &VtxSortedIndexVector )
{
    for (unsigned int i = 0; i < VtxSortedIndexVector.size(); i++) {
        int index = VtxSortedIndexVector[i];
        if (index < 0) continue;
        if ( (unsigned int) index == trueVtxIndex ) return i;
    }
    return -1;
}

std::vector<int>
ZMuMuValidationTreeMaker::isTightMuonWithoutVtxCut(const vector<edm::Ptr<pat::Muon>> &patMuons) 
{
    vector<int> muon_index;
    for ( unsigned int i = 0; i < patMuons.size(); i++) {
        Ptr<pat::Muon> pat_muon = patMuons[i];
        if(!pat_muon->isLooseMuon() || pat_muon->pt()<10. ) continue;
        if(!pat_muon->innerTrack().isNonnull()) continue;
        if(!pat_muon->globalTrack().isNonnull()) continue;
        if(pat_muon->globalTrack()->normalizedChi2() > 10. ) continue;
        if(pat_muon->globalTrack()->hitPattern().numberOfValidMuonHits() <= 0) continue;
        if(pat_muon->numberOfMatchedStations() <= 1) continue;
        if(pat_muon->innerTrack()->hitPattern().numberOfValidPixelHits() <= 0) continue;
        if(pat_muon->innerTrack()->hitPattern().trackerLayersWithMeasurement() <= 5) continue;
        if(pat_muon->isolationR03().sumPt/pat_muon->pt()>0.05) continue;
        muon_index.push_back(i);
    }

    return muon_index;
}

void
ZMuMuValidationTreeMaker::fillDescriptions( edm::ConfigurationDescriptions &descriptions )
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault( desc );
}


DEFINE_FWK_MODULE( ZMuMuValidationTreeMaker );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

