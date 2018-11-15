import FWCore.ParameterSet.Config as cms

def getDiPhotonSystematicsList():

    phosystlabels = []
    for direction in ["Up","Down"]:
        phosystlabels.append("MvaShift%s01sigma" % direction)
        phosystlabels.append("SigmaEOverEShift%s01sigma" % direction)
        phosystlabels.append("MaterialCentralBarrel%s01sigma" % direction)
        phosystlabels.append("MaterialOuterBarrel%s01sigma" % direction)
        phosystlabels.append("MaterialForward%s01sigma" % direction)
        phosystlabels.append("FNUFEB%s01sigma" % direction)
        phosystlabels.append("FNUFEE%s01sigma" % direction)
        phosystlabels.append("MCScaleGain6EB%s01sigma" % direction)
        phosystlabels.append("MCScaleGain1EB%s01sigma" % direction)
        for r9 in ["HighR9","LowR9"]:
            for region in ["EB","EE"]:
                phosystlabels.append("ShowerShape%s%s%s01sigma"%(r9,region,direction))
                phosystlabels.append("MCScale%s%s%s01sigma" % (r9,region,direction))
                for var in ["Rho","Phi"]:
                    phosystlabels.append("MCSmear%s%s%s%s01sigma" % (r9,region,var,direction))
        #variablesToUse.append("UnmatchedPUWeight%s01sigma[1,-999999.,999999.] := weight(\"UnmatchedPUWeight%s01sigma\")" % (direction,direction))
        #variablesToUse.append("MvaLinearSyst%s01sigma[1,-999999.,999999.] := weight(\"MvaLinearSyst%s01sigma\")" % (direction,direction))
        #variablesToUse.append("LooseMvaSF%s01sigma[1,-999999.,999999.] := weight(\"LooseMvaSF%s01sigma\")" % (direction,direction))
        #variablesToUse.append("PreselSF%s01sigma[1,-999999.,999999.] := weight(\"PreselSF%s01sigma\")" % (direction,direction))
        #variablesToUse.append("electronVetoSF%s01sigma[1,-999999.,999999.] := weight(\"electronVetoSF%s01sigma\")" % (direction,direction))
        #variablesToUse.append("TriggerWeight%s01sigma[1,-999999.,999999.] := weight(\"TriggerWeight%s01sigma\")" % (direction,direction))
        #variablesToUse.append("FracRVWeight%s01sigma[1,-999999.,999999.] := weight(\"FracRVWeight%s01sigma\")" % (direction,direction))
        #variablesToUse.append("FracRVNvtxWeight%s01sigma[1,-999999.,999999.] := weight(\"FracRVNvtxWeight%s01sigma\")" % (direction,direction))
        #variablesToUse.append("JetBTagCutWeight%s01sigma[1,-999999.,999999.] := weight(\"JetBTagCutWeight%s01sigma\")" % (direction,direction))
        #variablesToUse.append("JetBTagReshapeWeight%s01sigma[1,-999999.,999999.] := weight(\"JetBTagReshapeWeight%s01sigma\")" % (direction,direction))
    return phosystlabels

def customizeSystematicsForMC(process):
    photonSmearBins = getattr(process,'photonSmearBins',None)
    photonScaleUncertBins = getattr(process,'photonScaleUncertBins',None)
    for pset in process.flashggDiPhotonSystematics.SystMethods:
        if photonSmearBins and pset.Label.value().startswith("MCSmear"):
            pset.BinList = photonSmearBins
        elif photonScaleUncertBins and pset.Label.value().count("Scale"):
            pset.BinList = photonScaleUncertBins

def includeScale_Central_Systematics(process):
    customizeSystematicsForMC(process)

def includeScale_Central(process):
    # Keep default MC central value behavior, remove all up/down shifts
    customizeSystematicsForMC(process)
    vpsetlist = [process.flashggDiPhotonSystematics.SystMethods]
    vpsetlist += [process.flashggDiPhotonSystematics.SystMethods2D]
    for vpset in vpsetlist:
        for pset in vpset:
            if type(pset.NSigmas) == type(cms.vint32()):
                pset.NSigmas = cms.vint32() # Do not perform shift
            else:
                pset.NSigmas = cms.PSet( firstVar = cms.vint32(), secondVar = cms.vint32() ) # Do not perform shift - 2D case

def includeScale(process):
    # By default remove the systematic entirely (central value and shifts)
    # For scale: put in central value, but omit shifts
    # TODO: this is wrong for sigE/E and possibly others - check!

    photonScaleBinsData = getattr(process,'photonScaleBinsData',None)
    if hasattr(process,'photonScaleBinsData'):
        print photonScaleBinsData, process.photonScaleBinsData
    process.flashggDiPhotonSystematics.SystMethods = customizeVPSetForData(process.flashggDiPhotonSystematics.SystMethods, photonScaleBinsData)
    process.flashggDiPhotonSystematics.SystMethods2D = customizeVPSetForData(process.flashggDiPhotonSystematics.SystMethods2D, photonScaleBinsData)

def prepareflashggDiPhotonSystematicsTask(process, processType, doSystematics = False):

    from flashgg.Systematics.SystematicsCustomize import useEGMTools
    process.load("flashgg.Systematics.flashggDiPhotonSystematics_cfi")
    process.flashggPreselectedDiPhotons.src = cms.InputTag('flashggDiPhotonSystematics')

    SystTask = cms.Task(process.flashggDiPhotonSystematics)

    useEGMTools(process)

    if processType == 'sig':
        if doSystematics:
            includeScale_Central_Systematics(process)
        else:
            includeScale_Central(process)
    elif processType == 'bkg':
        includeScale_Central(process)
    elif processType == 'data':
        includeScale(process)
    else:
        print "Please choose processType which is 'sig', 'bkg', 'data' in prepareflashggDiPhotonSystematicsTask(..., processType, ...) "

    if doSystematics:
        for phosystlabel in getDiPhotonSystematicsList():
            setattr( process, 'flashggPreselectedDiPhotons' + phosystlabel,
                            process.flashggPreselectedDiPhotons.clone(
                                src = cms.InputTag('flashggDiPhotonSystematics', phosystlabel)
                            )
            )

            SystTask.add( getattr( process, 'flashggPreselectedDiPhotons' + phosystlabel ) )

    return SystTask
