
#   include "readTimeControls.H"

    // Remove unused variable warnings
    adjustTimeStep = adjustTimeStep;
    maxCo = maxCo;
    maxDeltaT = maxDeltaT;


#   include "readFluidPIMPLEControls.H"

    // Remove unused variable warnings
    methodUp = methodUp;

//     momentumPredictor = momentumPredictor;
//     correctPhi = correctPhi;
//     checkMeshCourantNo = checkMeshCourantNo;

    nOuterCorr = nOuterCorr;
    nCorr = nCorr;
    nNonOrthCorr = nNonOrthCorr;

    nMinOuterCorr = nMinOuterCorr;
    maxRelVolErr = maxRelVolErr;
    outerCorrBreak = outerCorrBreak;

    maxConvectiveCo = maxConvectiveCo;
    maxInterfaceCo = maxInterfaceCo;

    removeSwirl = removeSwirl;

    volCorrByHeight = volCorrByHeight;
    volCorrRelThreshold = volCorrRelThreshold;

    lorentzForceRotationalFactor = lorentzForceRotationalFactor;
    lorentzForceVolumeFactor = lorentzForceVolumeFactor;
