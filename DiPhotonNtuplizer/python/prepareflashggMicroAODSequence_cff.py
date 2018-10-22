import FWCore.ParameterSet.Config as cms

def prepareMicroAODTask (process, isData):
    process.load("myflashggPlugins.vbf_analysis.flashggDiphotonSequence_cff")
    if not isData:
        process.load("myflashggPlugins.vbf_analysis.flashggHTXS_cfi")
        process.load("flashgg.MicroAOD.flashggMicroAODGenSequence_cff")

    from myflashggPlugins.vbf_analysis.flashggJets_cfi import buildFinalJetsTask
    buildFinalJetsTask(process, isData)

    #setattr( process, 'MicroAODTask', cms.Task() )
    #getattr( process, 'MicroAODTask', cms.Task() ).add(*[getattr(process,prod) for prod in process.producers_()])
    #getattr( process, 'MicroAODTask', cms.Task() ).add(*[getattr(process,filt) for filt in process.filters_()])

    from PhysicsTools.SelectorUtils.tools.vid_id_tools import DataFormat,switchOnVIDPhotonIdProducer,setupAllVIDIdsInModule,setupVIDPhotonSelection
    switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD)
    #my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff']
    my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1_cff']
    for idmod in my_id_modules:
        setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)
    process.flashggPhotons.effAreasConfigFile = cms.FileInPath("RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_TrueVtx.txt")
    process.flashggPhotons.egmMvaValuesMap = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRunIIFall17v1Values")

    setattr( process, 'MicroAODTask', cms.Task() )
    getattr( process, 'MicroAODTask', cms.Task() ).add(*[getattr(process,prod) for prod in process.producers_()])
    getattr( process, 'MicroAODTask', cms.Task() ).add(*[getattr(process,filt) for filt in process.filters_()])

    #process.MicroAODTask.add(process.egmPhotonIDSequence)

    return process.MicroAODTask
