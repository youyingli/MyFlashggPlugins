#ifndef __JETSYSTEMATICS__
#define __JETSYSTEMATICS__

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "flashgg/DataFormats/interface/Jet.h"

#include <tuple>

namespace flashggAnalysisNtuplizer {

    const double 
    JECUncertainty ( const flashgg::Jet &jet, const edm::EventSetup &iSetup )
    {
        edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
        iSetup.get<JetCorrectionsRecord>().get("AK4PFchs", JetCorParColl); 
        JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
        JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);
    
        jecUnc->setJetEta(jet.eta());
        jecUnc->setJetPt (jet.pt()); 
        float unc = jecUnc->getUncertainty(true);
        delete jecUnc;
    
        return unc;
    }
    
    const std::tuple<float, float, float>
    JERUncertainty ( const flashgg::Jet &jet, double rho, const edm::EventSetup &iSetup )
    {
        JME::JetResolution resolution = JME::JetResolution::get(iSetup, "AK4PFchs_pt");
        JME::JetResolutionScaleFactor resolution_sf = JME::JetResolutionScaleFactor::get(iSetup, "AK4PFchs");
    
        JME::JetParameters params;
        params.setJetPt(jet.pt()).setJetEta(jet.eta()).setRho(rho);
    
        float res = resolution.getResolution(params);
        float sf      = resolution_sf.getScaleFactor(params);
        float sf_up   = resolution_sf.getScaleFactor(params, Variation::UP);
        float sf_down = resolution_sf.getScaleFactor(params, Variation::DOWN);
    
        float jet_pt = jet.pt();
        auto genjet = jet.genJet();
    
        if (   genjet != nullptr 
            && deltaR(jet, *(genjet)) < 0.2 
            && fabs(jet_pt - genjet->pt()) < 3 * res * jet_pt )
        {
            float genjet_pt = genjet->pt();
            float c_jer      = max<float>( 1 + (sf - 1)      * (jet_pt - genjet_pt) / jet_pt, 0. );
            float c_jer_up   = max<float>( 1 + (sf_up - 1)   * (jet_pt - genjet_pt) / jet_pt, 0. );
            float c_jer_down = max<float>( 1 + (sf_down - 1) * (jet_pt - genjet_pt) / jet_pt, 0. );
            return std::make_tuple( c_jer, c_jer_up, c_jer_down );
        } else {
            float rnd = jet.userFloat("rnd_g_JER");
            float c_jer      = ( 1. + rnd * std::sqrt(max<float>( sf * sf           - 1, 0.)) * res );
            float c_jer_up   = ( 1. + rnd * std::sqrt(max<float>( sf_up * sf_up     - 1, 0.)) * res );
            float c_jer_down = ( 1. + rnd * std::sqrt(max<float>( sf_down * sf_down - 1, 0.)) * res );
            return std::make_tuple( c_jer, c_jer_up, c_jer_down );
        }
    }
}

#endif
