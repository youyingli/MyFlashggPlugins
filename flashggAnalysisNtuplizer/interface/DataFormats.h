#ifndef __FLASHGG__DATAFORMATS__
#define __FLASHGG__DATAFORMATS__

#include <TTree.h>
#include <vector>

class flashggAnalysisTreeFormatStd
{
    public:

        int NPu                                       ;
        int NVtx                                      ;
        bool passTrigger                              ;
        float genweight                               ;
        float Rho                                     ;
        float PVz                                     ;
        float BSsigmaz                                ;
        bool Flag_HBHENoiseFilter                     ;
        bool Flag_HBHENoiseIsoFilter                  ;
        bool Flag_EcalDeadCellTriggerPrimitiveFilter  ;
        bool Flag_goodVertices                        ;
        bool Flag_globalSuperTightHalo2016Filter      ;
        bool Flag_BadPFMuonFilter                     ;
        bool Flag_BadChargedCandidateFilter           ;
        bool Flag_ecalBadCalibFilter                  ;
        bool Flag_eeBadScFilter                       ;
                                            
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
        bool  dipho_leadhasPixelSeed     ;
        bool  dipho_leadGenMatch         ;
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
        bool  dipho_subleadhasPixelSeed  ;
        bool  dipho_subleadGenMatch      ;
        int   dipho_subleadGenMatchType  ;
        float dipho_SelectedVz           ;
        float dipho_GenVz                ;
        
        int elecs_size;
        vector<int>   elecs_Charge              ; 
        vector<float> elecs_Pt                  ; 
        vector<float> elecs_Eta                 ; 
        vector<float> elecs_Phi                 ; 
        vector<float> elecs_Energy              ; 
        vector<float> elecs_EtaSC               ; 
        vector<float> elecs_PhiSC               ;
        vector<float> elecs_GsfTrackDz          ; 
        vector<float> elecs_GsfTrackDxy         ; 
        vector<bool>  elecs_EGMCutBasedIDVeto   ; 
        vector<bool>  elecs_EGMCutBasedIDLoose  ; 
        vector<bool>  elecs_EGMCutBasedIDMedium ; 
        vector<bool>  elecs_EGMCutBasedIDTight  ; 
        vector<bool>  elecs_fggPhoVeto          ;
        vector<bool>  elecs_tmpPhoVeto          ;
        vector<bool>  elecs_GenMatch            ;
        vector<int>   elecs_GenPdgID            ; 
        vector<float> elecs_GenPt               ;
        vector<float> elecs_GenEta              ;
        vector<float> elecs_GenPhi              ;

        int muons_size;
        vector<int>   muons_Charge                ; 
        vector<float> muons_MuonType              ; 
        vector<float> muons_Pt                    ; 
        vector<float> muons_Eta                   ; 
        vector<float> muons_Phi                   ; 
        vector<float> muons_Energy                ; 
        vector<float> muons_BestTrackDz           ; 
        vector<float> muons_BestTrackDxy          ; 
        vector<float> muons_PFIsoDeltaBetaCorrR04 ; 
        vector<float> muons_TrackerBasedIsoR03    ; 
        vector<bool>  muons_CutBasedIdMedium      ; 
        vector<bool>  muons_CutBasedIdTight       ;
        vector<bool>  muons_GenMatch              ;
        vector<int>   muons_GenPdgID              ;
        vector<float> muons_GenPt                 ;
        vector<float> muons_GenEta                ;
        vector<float> muons_GenPhi                ;

        int jets_size;
        vector<float> jets_Pt                                             ;
        vector<float> jets_Eta                                            ;
        vector<float> jets_Phi                                            ;
        vector<float> jets_Mass                                           ;
        vector<float> jets_Energy                                         ;
        vector<float> jets_PtRaw                                          ;
        vector<float> jets_QGL                                            ;
        vector<float> jets_RMS                                            ;
        vector<float> jets_puJetIdMVA                                     ;
        vector<bool>  jets_GenJetMatch                                    ;
        vector<float> jets_pfCombinedInclusiveSecondaryVertexV2BJetTags   ; 
        vector<float> jets_pfCombinedMVAV2BJetTags                        ; 
        vector<float> jets_pfDeepCSVJetTags_probb                         ; 
        vector<float> jets_pfDeepCSVJetTags_probbb                        ; 
        vector<float> jets_pfDeepCSVJetTags_probc                         ; 
        vector<float> jets_pfDeepCSVJetTags_probudsg                      ; 


        float met_Pt    ; 
        float met_Phi   ; 
        float met_Px    ; 
        float met_Py    ; 
        float met_SumET ; 

        int GenParticles_size;
        vector<float> GenParticles_Pt     ; 
        vector<float> GenParticles_Eta    ; 
        vector<float> GenParticles_Phi    ; 
        vector<float> GenParticles_Mass   ; 
        vector<int>   GenParticles_PdgID  ; 
        vector<int>   GenParticles_Status ; 
        vector<int>   GenParticles_nMo    ; 
        vector<int>   GenParticles_nDa    ; 

        int HTXSstage0cat ;
        int HTXSstage1cat ;
        int HTXSnjets     ;
        float HTXSpTH     ;
        float HTXSpTV     ;

        void Initialzation() {

            NPu                                     = -999;
            NVtx                                    = -999;
            passTrigger                             = false;
            genweight                               = -999.; 
            Rho                                     = -999.;
            PVz                                     = -999.;
            BSsigmaz                                = -999.;
            Flag_HBHENoiseFilter                    = false; 
            Flag_HBHENoiseIsoFilter                 = false; 
            Flag_EcalDeadCellTriggerPrimitiveFilter = false; 
            Flag_goodVertices                       = false; 
            Flag_globalSuperTightHalo2016Filter     = false; 
            Flag_BadPFMuonFilter                    = false; 
            Flag_BadChargedCandidateFilter          = false; 
            Flag_ecalBadCalibFilter                 = false; 
            Flag_eeBadScFilter                      = false; 

            dipho_mass                     = -999.;
            dipho_pt                       = -999.;
            dipho_leadPt                   = -999.; 
            dipho_leadEta                  = -999.; 
            dipho_leadPhi                  = -999.; 
            dipho_leadE                    = -999.; 
            dipho_leadsigEOverE            = -999.; 
            dipho_leadR9                   = -999.; 
            dipho_leadsieie                = -999.; 
            dipho_leadhoe                  = -999.; 
            dipho_leadIDMVA                = -999.;
            dipho_leadhasPixelSeed         = false;
            dipho_leadGenMatch             = false; 
            dipho_leadGenMatchType         = -999; 
            dipho_subleadPt                = -999.; 
            dipho_subleadEta               = -999.; 
            dipho_subleadPhi               = -999.; 
            dipho_subleadE                 = -999.; 
            dipho_subleadsigEOverE         = -999.; 
            dipho_subleadR9                = -999.; 
            dipho_subleadsieie             = -999.; 
            dipho_subleadhoe               = -999.; 
            dipho_subleadIDMVA             = -999.; 
            dipho_subleadhasPixelSeed      = false;
            dipho_subleadGenMatch          = false; 
            dipho_subleadGenMatchType      = -999;
            dipho_SelectedVz               = -999.;
            dipho_GenVz                    = -999.;

            elecs_size                  = -999;
            elecs_Charge                .clear(); 
            elecs_Pt                    .clear(); 
            elecs_Eta                   .clear(); 
            elecs_Phi                   .clear(); 
            elecs_Energy                .clear(); 
            elecs_EtaSC                 .clear(); 
            elecs_PhiSC                 .clear();
            elecs_GsfTrackDz            .clear();
            elecs_GsfTrackDxy           .clear();
            elecs_EGMCutBasedIDVeto     .clear(); 
            elecs_EGMCutBasedIDLoose    .clear(); 
            elecs_EGMCutBasedIDMedium   .clear(); 
            elecs_EGMCutBasedIDTight    .clear(); 
            elecs_fggPhoVeto            .clear();
            elecs_tmpPhoVeto            .clear();
            elecs_GenMatch              .clear(); 
            elecs_GenPdgID              .clear(); 
            elecs_GenPt                 .clear(); 
            elecs_GenEta                .clear(); 
            elecs_GenPhi                .clear(); 

            muons_size                  = -999;
            muons_Charge                .clear(); 
            muons_MuonType              .clear(); 
            muons_Pt                    .clear(); 
            muons_Eta                   .clear(); 
            muons_Phi                   .clear(); 
            muons_Energy                .clear(); 
            muons_BestTrackDz           .clear(); 
            muons_BestTrackDxy          .clear(); 
            muons_PFIsoDeltaBetaCorrR04 .clear(); 
            muons_TrackerBasedIsoR03    .clear(); 
            muons_CutBasedIdMedium      .clear(); 
            muons_CutBasedIdTight       .clear();
            muons_GenMatch              .clear(); 
            muons_GenPdgID              .clear(); 
            muons_GenPt                 .clear(); 
            muons_GenEta                .clear(); 
            muons_GenPhi                .clear(); 
 
            jets_size   = -999;
            jets_Pt                                           .clear(); 
            jets_Eta                                          .clear(); 
            jets_Phi                                          .clear(); 
            jets_Mass                                         .clear(); 
            jets_Energy                                       .clear(); 
            jets_PtRaw                                        .clear(); 
            jets_QGL                                          .clear(); 
            jets_RMS                                          .clear(); 
            jets_puJetIdMVA                                   .clear(); 
            jets_GenJetMatch                                  .clear();
            jets_pfCombinedInclusiveSecondaryVertexV2BJetTags .clear(); 
            jets_pfCombinedMVAV2BJetTags                      .clear(); 
            jets_pfDeepCSVJetTags_probb                       .clear(); 
            jets_pfDeepCSVJetTags_probbb                      .clear(); 
            jets_pfDeepCSVJetTags_probc                       .clear(); 
            jets_pfDeepCSVJetTags_probudsg                    .clear(); 

            met_Pt    = -999.;
            met_Phi   = -999.; 
            met_Px    = -999.; 
            met_Py    = -999.;
            met_SumET = -999.;

            GenParticles_size = -999;
            GenParticles_Pt     .clear(); 
            GenParticles_Eta    .clear(); 
            GenParticles_Phi    .clear(); 
            GenParticles_Mass   .clear(); 
            GenParticles_PdgID  .clear(); 
            GenParticles_Status .clear(); 
            GenParticles_nMo    .clear(); 
            GenParticles_nDa    .clear(); 

            HTXSstage0cat = -999; 
            HTXSstage1cat = -999; 
            HTXSnjets     = -999; 
            HTXSpTH       = -999.; 
            HTXSpTV       = -999.; 
        }

        void RegisterTree (TTree* Tree) {

            tree = Tree;

            tree->Branch( "EvtInfo.NPu"                                     , &NPu                                     , "EvtInfo.NPu/I"            );
            tree->Branch( "EvtInfo.NVtx"                                    , &NVtx                                    , "EvtInfo.NVtx/I"           );
            tree->Branch( "EvtInfo.passTrigger"                             , &passTrigger                             , "EvtInfo.passTrigger/O"    );
            tree->Branch( "EvtInfo.genweight"                               , &genweight                               , "EvtInfo.genweight/F"      );
            tree->Branch( "EvtInfo.Rho"                                     , &Rho                                     , "EvtInfo.Rho/F"            );
            tree->Branch( "EvtInfo.PVz"                                     , &PVz                                     , "EvtInfo.PVz/F"            );
            tree->Branch( "EvtInfo.BSsigmaz"                                , &BSsigmaz                                , "EvtInfo.BSsigmaz/F"       );
            tree->Branch( "EvtInfo.Flag_HBHENoiseFilter"                    , &Flag_HBHENoiseFilter                    , "EvtInfo.Flag_HBHENoiseFilter/O"                         );  
            tree->Branch( "EvtInfo.Flag_HBHENoiseIsoFilter"                 , &Flag_HBHENoiseIsoFilter                 , "EvtInfo.Flag_HBHENoiseIsoFilter/O"                      );  
            tree->Branch( "EvtInfo.Flag_EcalDeadCellTriggerPrimitiveFilter" , &Flag_EcalDeadCellTriggerPrimitiveFilter , "EvtInfo.Flag_EcalDeadCellTriggerPrimitiveFilter/O"      );  
            tree->Branch( "EvtInfo.Flag_goodVertices"                       , &Flag_goodVertices                       , "EvtInfo.Flag_goodVertices/O"                            );  
            tree->Branch( "EvtInfo.Flag_globalSuperTightHalo2016Filter"     , &Flag_globalSuperTightHalo2016Filter     , "EvtInfo.Flag_globalSuperTightHalo2016Filter/O"          );  
            tree->Branch( "EvtInfo.Flag_BadPFMuonFilter"                    , &Flag_BadPFMuonFilter                    , "EvtInfo.Flag_BadPFMuonFilter/O"                         );  
            tree->Branch( "EvtInfo.Flag_BadChargedCandidateFilter"          , &Flag_BadChargedCandidateFilter          , "EvtInfo.Flag_BadChargedCandidateFilter/O"               );  
            tree->Branch( "EvtInfo.Flag_ecalBadCalibFilter"                 , &Flag_ecalBadCalibFilter                 , "EvtInfo.Flag_ecalBadCalibFilter/O"                      );  
            tree->Branch( "EvtInfo.Flag_eeBadScFilter"                      , &Flag_eeBadScFilter                      , "EvtInfo.Flag_eeBadScFilter/O"                           );  

            tree->Branch( "DiPhoInfo.mass"                    , &dipho_mass                   , "DiPhoInfo.mass/F"                 );
            tree->Branch( "DiPhoInfo.pt"                      , &dipho_pt                     , "DiPhoInfo.pt/F"                   );
            tree->Branch( "DiPhoInfo.leadPt"                  , &dipho_leadPt                 , "DiPhoInfo.leadPt/F"               );
            tree->Branch( "DiPhoInfo.leadEta"                 , &dipho_leadEta                , "DiPhoInfo.leadEta/F"              );
            tree->Branch( "DiPhoInfo.leadPhi"                 , &dipho_leadPhi                , "DiPhoInfo.leadPhi/F"              );
            tree->Branch( "DiPhoInfo.leadE"                   , &dipho_leadE                  , "DiPhoInfo.leadE/F"                );
            tree->Branch( "DiPhoInfo.leadsigEOverE"           , &dipho_leadsigEOverE          , "DiPhoInfo.leadsigEOverE/F"        );
            tree->Branch( "DiPhoInfo.leadR9"                  , &dipho_leadR9                 , "DiPhoInfo.leadR9/F"               );
            tree->Branch( "DiPhoInfo.leadsieie"               , &dipho_leadsieie              , "DiPhoInfo.leadsieie/F"            );
            tree->Branch( "DiPhoInfo.leadhoe"                 , &dipho_leadhoe                , "DiPhoInfo.leadhoe/F"              );
            tree->Branch( "DiPhoInfo.leadIDMVA"               , &dipho_leadIDMVA              , "DiPhoInfo.leadIDMVA/F"            );
            tree->Branch( "DiPhoInfo.leadhasPixelSeed"        , &dipho_leadhasPixelSeed       , "DiPhoInfo.leadhasPixelSeed/O"     );
            tree->Branch( "DiPhoInfo.leadGenMatch"            , &dipho_leadGenMatch           , "DiPhoInfo.leadGenMatch/O"         );
            tree->Branch( "DiPhoInfo.leadGenMatchType"        , &dipho_leadGenMatchType       , "DiPhoInfo.leadGenMatchType/I"     );
            tree->Branch( "DiPhoInfo.subleadPt"               , &dipho_subleadPt              , "DiPhoInfo.subleadPt/F"            );
            tree->Branch( "DiPhoInfo.subleadEta"              , &dipho_subleadEta             , "DiPhoInfo.subleadEta/F"           );
            tree->Branch( "DiPhoInfo.subleadPhi"              , &dipho_subleadPhi             , "DiPhoInfo.subleadPhi/F"           );
            tree->Branch( "DiPhoInfo.subleadE"                , &dipho_subleadE               , "DiPhoInfo.subleadE/F"             );
            tree->Branch( "DiPhoInfo.subleadsigEOverE"        , &dipho_subleadsigEOverE       , "DiPhoInfo.subleadsigEOverE/F"     );
            tree->Branch( "DiPhoInfo.subleadR9"               , &dipho_subleadR9              , "DiPhoInfo.subleadR9/F"            );
            tree->Branch( "DiPhoInfo.subleadsieie"            , &dipho_subleadsieie           , "DiPhoInfo.subleadsieie/F"         );
            tree->Branch( "DiPhoInfo.subleadhoe"              , &dipho_subleadhoe             , "DiPhoInfo.subleadhoe/F"           );
            tree->Branch( "DiPhoInfo.subleadIDMVA"            , &dipho_subleadIDMVA           , "DiPhoInfo.subleadIDMVA/F"         );
            tree->Branch( "DiPhoInfo.subleadhasPixelSeed"     , &dipho_subleadhasPixelSeed    , "DiPhoInfo.subleadhasPixelSeed/O"  );
            tree->Branch( "DiPhoInfo.subleadGenMatch"         , &dipho_subleadGenMatch        , "DiPhoInfo.subleadGenMatch/O"      );
            tree->Branch( "DiPhoInfo.subleadGenMatchType"     , &dipho_subleadGenMatchType    , "DiPhoInfo.subleadGenMatchType/I"  );
            tree->Branch( "DiPhoInfo.SelectedVz"              , &dipho_SelectedVz             , "DiPhoInfo.SelectedVz/F"           );
            tree->Branch( "DiPhoInfo.GenVz"                   , &dipho_GenVz                  , "DiPhoInfo.GenVz/F"                );

            tree->Branch( "ElecInfo.Size"                        , &elecs_size                     , "ElecInfo.Size/I" );
            tree->Branch( "ElecInfo.Charge"                      , &elecs_Charge                ); 
            tree->Branch( "ElecInfo.Pt"                          , &elecs_Pt                    ); 
            tree->Branch( "ElecInfo.Eta"                         , &elecs_Eta                   ); 
            tree->Branch( "ElecInfo.Phi"                         , &elecs_Phi                   ); 
            tree->Branch( "ElecInfo.Energy"                      , &elecs_Energy                ); 
            tree->Branch( "ElecInfo.EtaSC"                       , &elecs_EtaSC                 ); 
            tree->Branch( "ElecInfo.PhiSC"                       , &elecs_PhiSC                 ); 
            tree->Branch( "ElecInfo.GsfTrackDz"                  , &elecs_GsfTrackDz            ); 
            tree->Branch( "ElecInfo.GsfTrackDxy"                 , &elecs_GsfTrackDxy           ); 
            tree->Branch( "ElecInfo.EGMCutBasedIDVeto"           , &elecs_EGMCutBasedIDVeto     ); 
            tree->Branch( "ElecInfo.EGMCutBasedIDLoose"          , &elecs_EGMCutBasedIDLoose    ); 
            tree->Branch( "ElecInfo.EGMCutBasedIDMedium"         , &elecs_EGMCutBasedIDMedium   ); 
            tree->Branch( "ElecInfo.EGMCutBasedIDTight"          , &elecs_EGMCutBasedIDTight    ); 
            tree->Branch( "ElecInfo.fggPhoVeto"                  , &elecs_fggPhoVeto            );
            tree->Branch( "ElecInfo.tmpPhoVeto"                  , &elecs_tmpPhoVeto            );
            tree->Branch( "ElecInfo.GenMatch"                    , &elecs_GenMatch              ); 
            tree->Branch( "ElecInfo.GenPdgID"                    , &elecs_GenPdgID              ); 
            tree->Branch( "ElecInfo.GenPt"                       , &elecs_GenPt                 ); 
            tree->Branch( "ElecInfo.GenEta"                      , &elecs_GenEta                ); 
            tree->Branch( "ElecInfo.GenPhi"                      , &elecs_GenPhi                ); 

            tree->Branch( "MuonInfo.Size"                        , &muons_size                     , "MuonInfo.Size/I" );
            tree->Branch( "MuonInfo.Charge"                      , &muons_Charge                ); 
            tree->Branch( "MuonInfo.MuonType"                    , &muons_MuonType              ); 
            tree->Branch( "MuonInfo.Pt"                          , &muons_Pt                    ); 
            tree->Branch( "MuonInfo.Eta"                         , &muons_Eta                   ); 
            tree->Branch( "MuonInfo.Phi"                         , &muons_Phi                   ); 
            tree->Branch( "MuonInfo.Energy"                      , &muons_Energy                ); 
            tree->Branch( "MuonInfo.BestTrackDz"                 , &muons_BestTrackDz           ); 
            tree->Branch( "MuonInfo.BestTrackDxy"                , &muons_BestTrackDxy          ); 
            tree->Branch( "MuonInfo.PFIsoDeltaBetaCorrR04"       , &muons_PFIsoDeltaBetaCorrR04 ); 
            tree->Branch( "MuonInfo.TrackerBasedIsoR03"          , &muons_TrackerBasedIsoR03    ); 
            tree->Branch( "MuonInfo.CutBasedIdMedium"            , &muons_CutBasedIdMedium      ); 
            tree->Branch( "MuonInfo.CutBasedIdTight"             , &muons_CutBasedIdTight       );
            tree->Branch( "MuonInfo.GenMatch"                    , &muons_GenMatch              ); 
            tree->Branch( "MuonInfo.GenPdgID"                    , &muons_GenPdgID              ); 
            tree->Branch( "MuonInfo.GenPt"                       , &muons_GenPt                 ); 
            tree->Branch( "MuonInfo.GenEta"                      , &muons_GenEta                ); 
            tree->Branch( "MuonInfo.GenPhi"                      , &muons_GenPhi                ); 

            tree->Branch( "jets_size"                                             , &jets_size      , "jets_size/I" );
            tree->Branch( "JetInfo.Pt"                                            , &jets_Pt                                           );
            tree->Branch( "JetInfo.Eta"                                           , &jets_Eta                                          );
            tree->Branch( "JetInfo.Phi"                                           , &jets_Phi                                          );
            tree->Branch( "JetInfo.Mass"                                          , &jets_Mass                                         );
            tree->Branch( "JetInfo.Energy"                                        , &jets_Energy                                       );
            tree->Branch( "JetInfo.PtRaw"                                         , &jets_PtRaw                                        );
            tree->Branch( "JetInfo.QGL"                                           , &jets_QGL                                          );
            tree->Branch( "JetInfo.RMS"                                           , &jets_RMS                                          );
            tree->Branch( "JetInfo.puJetIdMVA"                                    , &jets_puJetIdMVA                                   );
            tree->Branch( "JetInfo.GenJetMatch"                                   , &jets_GenJetMatch                                  );
            tree->Branch( "JetInfo.pfCombinedInclusiveSecondaryVertexV2BJetTags"  , &jets_pfCombinedInclusiveSecondaryVertexV2BJetTags ); 
            tree->Branch( "JetInfo.pfCombinedMVAV2BJetTags"                       , &jets_pfCombinedMVAV2BJetTags                      ); 
            tree->Branch( "JetInfo.pfDeepCSVJetTags_probb"                        , &jets_pfDeepCSVJetTags_probb                       ); 
            tree->Branch( "JetInfo.pfDeepCSVJetTags_probbb"                       , &jets_pfDeepCSVJetTags_probbb                      ); 
            tree->Branch( "JetInfo.pfDeepCSVJetTags_probc"                        , &jets_pfDeepCSVJetTags_probc                       ); 
            tree->Branch( "JetInfo.pfDeepCSVJetTags_probudsg"                     , &jets_pfDeepCSVJetTags_probudsg                    ); 

            tree->Branch( "met_Pt"          , &met_Pt       , "met_Pt/F"    );  
            tree->Branch( "met_Phi"         , &met_Phi      , "met_Phi/F"   ); 
            tree->Branch( "met_Px"          , &met_Px       , "met_Px/F"    ); 
            tree->Branch( "met_Py"          , &met_Py       , "met_Py/F"    ); 
            tree->Branch( "met_SumET"       , &met_SumET    , "met_SumET/F" );

            tree->Branch( "GenPartInfo.size"   , &GenParticles_size , "GenPartInfo.size/I" );
            tree->Branch( "GenPartInfo.Pt"     , &GenParticles_Pt     ); 
            tree->Branch( "GenPartInfo.Eta"    , &GenParticles_Eta    ); 
            tree->Branch( "GenPartInfo.Phi"    , &GenParticles_Phi    ); 
            tree->Branch( "GenPartInfo.Mass"   , &GenParticles_Mass   ); 
            tree->Branch( "GenPartInfo.PdgID"  , &GenParticles_PdgID  ); 
            tree->Branch( "GenPartInfo.Status" , &GenParticles_Status ); 
            tree->Branch( "GenPartInfo.nMo"    , &GenParticles_nMo    ); 
            tree->Branch( "GenPartInfo.nDa"    , &GenParticles_nDa    ); 

            tree->Branch( "HTXSstage0cat"   , &HTXSstage0cat  , "HTXSstage0cat/I" ); 
            tree->Branch( "HTXSstage1cat"   , &HTXSstage1cat  , "HTXSstage1cat/I" ); 
            tree->Branch( "HTXSnjets"       , &HTXSnjets      , "HTXSnjets/I"     ); 
            tree->Branch( "HTXSpTH"         , &HTXSpTH        , "HTXSpTH/F"       ); 
            tree->Branch( "HTXSpTV"         , &HTXSpTV        , "HTXSpTV/F"       ); 

        }

        void InfoFill() { tree->Fill(); }
    private:
        TTree* tree;

};

#endif