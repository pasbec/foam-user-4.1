
#   include "readTimeControls.H"

    // Remove unused variable warnings
    adjustTimeStep = adjustTimeStep;
    maxCo = maxCo;
    maxDeltaT = maxDeltaT;


#   include "readPIMPLEControls.H"

    // Remove unused variable warnings
    nOuterCorr = nOuterCorr;
    nCorr = nCorr;
    nNonOrthCorr = nNonOrthCorr;
    momentumPredictor = momentumPredictor;
    transonic = transonic;


#   include "readExtendedPIMPLEControls.H"

    // Remove unused variable warnings
    nMinOuterCorr = nMinOuterCorr;
    maxRelVolErr = maxRelVolErr;
    maxInterfaceCo = maxInterfaceCo;
    maxConvectiveCo = maxConvectiveCo;
    outerCorrBreak = outerCorrBreak;
    removeSwirl = removeSwirl;
    correctVolByHeight = correctVolByHeight;


#   include "readLorentzForceControls.H"

    // Remove unused variable warnings
    lorentzForceRotationalFactor = lorentzForceRotationalFactor;
    lorentzForceVolumeFactor = lorentzForceVolumeFactor;
    lorentzForceNonAlberto = lorentzForceNonAlberto;
    lorentzForceUpdateSwitch = lorentzForceUpdateSwitch;
    lorentzForceUpdateOutputCycle = lorentzForceUpdateOutputCycle;
    lorentzForceUpdateTimeStepCycle = lorentzForceUpdateTimeStepCycle;
    lorentzForceUpdateTimeCycle = lorentzForceUpdateTimeCycle;
    lorentzForceUpdateUtility = lorentzForceUpdateUtility;