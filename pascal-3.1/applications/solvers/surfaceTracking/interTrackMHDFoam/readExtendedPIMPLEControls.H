
    int nMinOuterCorr =
        pimple.lookupOrDefault<int>("nMinOuterCorr", 10);

    scalar maxRelVolErr =
        pimple.lookupOrDefault<scalar>("maxRelVolErr", 1e-08);

    scalar maxInterfaceCo =
        pimple.lookupOrDefault<scalar>("maxInterfaceCo", 1.0);
    scalar maxConvectiveCo = readScalar(runTime.controlDict().lookup("maxCo"));

    bool outerCorrBreak = false;

    label removeSwirl =
          pimple.lookupOrDefault<label>("removeSwirl", -1);

    label correctVolByHeight =
          pimple.lookupOrDefault<label>("correctVolByHeight", -1);

//     Switch correctPhi =
//          pimple.lookupOrDefault<bool>("correctPhi", false);

//     Switch checkMeshCourantNo =
//          pimple.lookupOrDefault<bool>("checkMeshCourantNo", true);