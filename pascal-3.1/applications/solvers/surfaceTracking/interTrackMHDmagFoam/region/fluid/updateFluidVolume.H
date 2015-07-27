
// TODO [High]: What about boundary points? If mass loss is too large, boundary points might mis-align...

    scalar curFluidVolumePhase1 = fvc::domainIntegrate(fluidIndicator).value();

    scalar curRelFluidVolumePhase1 = curFluidVolumePhase1/(iniFluidVolumePhase1+VSMALL);

    fluidVolumeCurrentSubDict.set<scalar>("phase1", curFluidVolumePhase1);

    Info << "Absolute volume of phase1: " << curFluidVolumePhase1 << " ( " << (100*curRelFluidVolumePhase1) << " % )" << nl;

// TODO [Low]: Volume correction with twoFluids?

//     if (interface.twoFluids())
//     {
//         scalar curFluidVolumePhase2 = fvc::domainIntegrate(1.0-fluidIndicator).value();
//         scalar curFluidVolumeTotal = curFluidVolumePhase1 + curFluidVolumePhase2;
//
//         fluidVolumeCurrentSubDict.set<scalar>("phase2", curFluidVolumePhase2);
//         fluidVolumeCurrentSubDict.set<scalar>("total", curFluidVolumeTotal);
//
//         Info << "Absolute volume of phase2: " << curFluidVolumePhase2 << " ( " << (100*curFluidVolumePhase2/(iniFluidVolumePhase2+VSMALL)) << " % )" << nl;
//         Info << "Absolute volume total: " << curFluidVolumeTotal << " ( " << (100*curFluidVolumeTotal/(iniFluidVolumeTotal+VSMALL)) << " % )" << nl;
//     }

    if
    (
        volCorrByHeight > -1
     && abs(curRelFluidVolumePhase1-1.0) >= volCorrRelThreshold
    )
    {
        label aPatchID = interface.aPatchID();

        scalar volCorr = curFluidVolumePhase1/iniFluidVolumePhase1;

        scalarField volCorrDeltaH = mesh.C().boundaryField()[aPatchID].component(volCorrByHeight)
                                    * (1.0 - volCorr);

        scalarField volCorrPhi = mesh.magSf().boundaryField()[aPatchID]
                                * volCorrDeltaH
                                / runTime.deltaT().value()
                                + fvc::meshPhi(U)().boundaryField()[aPatchID];

        interface.movePoints(volCorrPhi);

        curFluidVolumePhase1 = fvc::domainIntegrate(fluidIndicator).value();

        curRelFluidVolumePhase1 = curFluidVolumePhase1/(iniFluidVolumePhase1+VSMALL);

        fluidVolumeCurrentSubDict.set<scalar>("phase1", curFluidVolumePhase1);

        Info << "Corrected volume of phase1!" << endl;

        Info << "New absolute volume of phase1: "
             << curFluidVolumePhase1 << " ( "
             << (100*curRelFluidVolumePhase1) << " % )" << nl;
    }

    Info << endl;
