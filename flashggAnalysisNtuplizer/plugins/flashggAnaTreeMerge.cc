#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "TTree.h"

#include "MyFlashggPlugins/flashggAnalysisNtuplizer/interface/flashggAnaTreeMakerWithSyst.h"

using namespace std;
using namespace edm;

class flashggAnaTreeMerge : public edm::EDAnalyzer
{
    public:
        explicit flashggAnaTreeMerge( const edm::ParameterSet & );
        ~flashggAnaTreeMerge();
    
        static void fillDescriptions( edm::ConfigurationDescriptions &descriptions );
    
    
    private:
    
        edm::Service<TFileService> fs_;
    
        virtual void beginJob() override;
        virtual void analyze( const edm::Event &, const edm::EventSetup & ) override;
        virtual void endJob() override;
    
    
        std::vector<std::string> DiphoSystNames_; 
        std::vector<flashggAnaTreeMakerWithSyst*> TreeMakeList_;

};

flashggAnaTreeMerge::flashggAnaTreeMerge( const edm::ParameterSet &iConfig ):
    DiphoSystNames_(iConfig.getParameter<vector<string>>( "diphosystnames" ))
{
    std::vector<edm::InputTag> diphotons        = iConfig.getParameter<std::vector<edm::InputTag>>( "diphotons" );
    std::vector<edm::InputTag> diphotonMVAs     = iConfig.getParameter<std::vector<edm::InputTag>>( "diphotonMVAs" );
    const auto&           NonDiphoSetting       = iConfig.getParameter<edm::ParameterSet>( "nondiphosetting" );

    for ( unsigned int i_diphoton = 0; i_diphoton < diphotons.size(); i_diphoton++ )
        TreeMakeList_.emplace_back( new flashggAnaTreeMakerWithSyst( diphotons[i_diphoton], diphotonMVAs[i_diphoton], 
                                                                     NonDiphoSetting, consumesCollector() ) );
}

flashggAnaTreeMerge::~flashggAnaTreeMerge()
{
    for ( const auto& TreeMake : TreeMakeList_ ) {
        delete TreeMake;
    }
}

void
flashggAnaTreeMerge::beginJob()
{
    int i_syst = 0;
    for ( auto& TreeMake : TreeMakeList_ ) {
        TTree* tree = fs_->make<TTree>( Form("flashggStdTree%s",DiphoSystNames_[i_syst].c_str()), "" );
        TreeMake->RegisterTree( tree );
        i_syst++;
    }
}

void
flashggAnaTreeMerge::endJob()
{
}

void
flashggAnaTreeMerge::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup )
{
    int itree = 0;
    for ( const auto& TreeMake : TreeMakeList_ ) {
        if (itree == 0) TreeMake->Analyze(iEvent, iSetup, false);
        else TreeMake->Analyze(iEvent, iSetup, true);
        itree++;
    }

}

void
flashggAnaTreeMerge::fillDescriptions( edm::ConfigurationDescriptions &descriptions )
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault( desc );
}

DEFINE_FWK_MODULE( flashggAnaTreeMerge );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
