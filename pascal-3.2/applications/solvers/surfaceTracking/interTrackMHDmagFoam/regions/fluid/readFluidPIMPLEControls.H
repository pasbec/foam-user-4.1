
    dictionary pimple =
        mesh_[fluidRegionID].solutionDict().subDict("PIMPLE");

    word methodUp =
        pimple.lookupOrDefault<word>("methodUp", "default");
    if
    (
       !(
            methodUp == "default"
         || methodUp == "alberto"
         || methodUp == "coupled"
        )
    )
    {
        Info << "Pressure veolcity Method 'methodUp' = " << methodUp
            << " is not supported!  Valid methods are:" << endl
            << "3" << endl
            << "(" << endl
            << "    default" << endl
            << "    alberto" << endl
            << "    coupled" << endl
            << ")"
            << abort(FatalError);
    }

//     bool momentumPredictor =
//         pimple.lookupOrDefault<bool>("momentumPredictor", true);
//
//     Switch correctPhi =
//          pimple.lookupOrDefault<bool>("correctPhi", false);
//
//     Switch checkMeshCourantNo =
//          pimple.lookupOrDefault<bool>("checkMeshCourantNo", true);


    // Basic PIMPLE loop settings

        int nOuterCorr(readInt(pimple.lookup("nOuterCorrectors")));

        int nCorr(readInt(pimple.lookup("nCorrectors")));

        int nNonOrthCorr =
            pimple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);


    // Pimple loop break critera

        int nMinOuterCorr =
            pimple.lookupOrDefault<label>("nMinOuterCorr", 10);

        bool outerCorrBreak = false;

        scalar maxRelVolErr =
            pimple.lookupOrDefault<scalar>("maxRelVolErr", 1e-08);


    // Courant numbers

        scalar maxConvectiveCo =
            pimple.lookupOrDefault<scalar>
            (
                "maxInterfaceCo",
                readScalar(runTime.controlDict().lookup("maxCo"))
            );
        scalar maxInterfaceCo =
            pimple.lookupOrDefault<scalar>("maxInterfaceCo", 1.0);


    // Remove swirl or one volocity component

        int removeSwirl =
            pimple.lookupOrDefault<label>("removeSwirl", -1);


    // Volume correction
    dictionary& volCorrSolutionDict(pimple.subDict("volumeCorrectionByHeight"));

        int volCorrByHeight =
            volCorrSolutionDict.lookupOrDefault<label>("heightCoord", -1);

        scalar volCorrRelThreshold =
            volCorrSolutionDict.lookupOrDefault<scalar>("relThreshold", 0.0);


    // Lorentz force settings
    dictionary& lorentzForceSolutionDict(pimple.subDict("lorentzForce"));

        scalar lorentzForceRotationalFactor =
            lorentzForceSolutionDict.lookupOrDefault<scalar>("rotationalFactor", 0.0);

        scalar lorentzForceVolumeFactor =
            lorentzForceSolutionDict.lookupOrDefault<scalar>("volumeFactor", 1.0);
