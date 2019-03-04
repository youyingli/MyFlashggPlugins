#ifndef __FLASHGGANATREEMAKERWITHSYST__
#define __FLASHGGANATREEMAKERWITHSYST__

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/DiPhotonMVAResult.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/DataFormats/interface/Met.h"
#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "MyFlashggPlugins/flashggAnalysisNtuplizer/interface/DataFormats.h"
#include "TTree.h"


class flashggAnaTreeMakerWithSyst
{
    public:
        flashggAnaTreeMakerWithSyst( const edm::InputTag &, const edm::InputTag &, const edm::ParameterSet &, edm::ConsumesCollector&& );
        ~flashggAnaTreeMakerWithSyst();
        void Analyze( const edm::Event &, const edm::EventSetup &, bool );
        void RegisterTree( TTree* );
    
    private:
    
        flashggAnalysisTreeFormatStd dataformat;
    
        edm::EDGetTokenT<edm::View<flashgg::DiPhotonCandidate>> diphotonToken_;
        edm::EDGetTokenT<edm::View<flashgg::DiPhotonMVAResult>> diphotonMVAToken_;
        std::vector<edm::InputTag> inputTagJets_;
        std::vector<edm::EDGetTokenT<edm::View<flashgg::Jet>>>  tokenJets_;
        edm::EDGetTokenT<edm::View<flashgg::Electron>>          electronToken_;
        edm::EDGetTokenT<edm::View<flashgg::Muon>>              muonToken_;
        edm::EDGetTokenT<edm::View<flashgg::Met>>               metToken_;
        edm::EDGetTokenT<edm::View<reco::Vertex>>               vertexToken_;
        edm::EDGetTokenT<reco::BeamSpot>                        beamSpotToken_;
        edm::EDGetTokenT<double>                                rhoTaken_;
        edm::EDGetTokenT<edm::View<reco::GenParticle>>          genParticleToken_;
        edm::EDGetTokenT<GenEventInfoProduct>                   genEventInfoToken_;
        edm::EDGetTokenT<edm::View<PileupSummaryInfo>>          pileUpToken_;
        edm::EDGetTokenT<edm::TriggerResults>                   triggerToken_;
        edm::EDGetTokenT<edm::TriggerResults>                   mettriggerToken_;
        edm::EDGetTokenT<HTXS::HiggsClassification>             newHTXSToken_;
        std::string pathName_;
        bool isMiniAOD_;
        bool storeSyst_;
        bool doHTXS_;
    
        typedef std::vector<edm::Handle<edm::View<flashgg::Jet>>> JetCollectionVector;
        edm::Handle<edm::View<flashgg::DiPhotonCandidate> >   diphotons;
        edm::Handle<edm::View<flashgg::DiPhotonMVAResult> >   diphotonMVAs;
        edm::Handle<edm::View<flashgg::Electron> >            electrons;
        edm::Handle<edm::View<flashgg::Muon> >                muons;
        edm::Handle<edm::View<flashgg::Met> >                 met;
        edm::Handle<edm::View<reco::Vertex> >                 primaryVertices;
        edm::Handle<reco::BeamSpot >                          recoBeamSpotHandle;
        edm::Handle<double >                                  rho;
        edm::Handle<edm::View<reco::GenParticle> >            genParticles;
        edm::Handle<GenEventInfoProduct >                     genEventInfo;
        edm::Handle<edm::View< PileupSummaryInfo> >           pileupInfo;
        edm::Handle<edm::TriggerResults >                     triggerHandle;
        edm::Handle<edm::TriggerResults >                     mettriggerHandle;
        edm::Handle<HTXS::HiggsClassification >               htxsClassification;
};

#endif
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
