import FWCore.ParameterSet.Config as cms

def myflashggCustomize(process, runMiniAOD):
    if runMiniAOD:
        process.flashggDiPhotons.useZerothVertexFromMicro = cms.bool(True)
