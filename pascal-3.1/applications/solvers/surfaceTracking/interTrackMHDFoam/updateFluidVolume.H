
    scalar curFluidVolumePhase1 = fvc::domainIntegrate(fluidIndicator).value();

    fluidVolumeCurrentSubDict.set<scalar>("phase1", curFluidVolumePhase1);

    Info << "Absolute volume of phase1: " << curFluidVolumePhase1 << " ( " << (100*curFluidVolumePhase1/(iniFluidVolumePhase1+VSMALL)) << " % )" << nl;

    if (interface.twoFluids())
    {
        scalar curFluidVolumePhase2 = fvc::domainIntegrate(1.0-fluidIndicator).value();
        scalar curFluidVolumeTotal = curFluidVolumePhase1 + curFluidVolumePhase2;

        fluidVolumeCurrentSubDict.set<scalar>("phase2", curFluidVolumePhase2);
        fluidVolumeCurrentSubDict.set<scalar>("total", curFluidVolumeTotal);

        Info << "Absolute volume of phase2: " << curFluidVolumePhase2 << " ( " << (100*curFluidVolumePhase2/(iniFluidVolumePhase2+VSMALL)) << " % )" << nl;
        Info << "Absolute volume total: " << curFluidVolumeTotal << " ( " << (100*curFluidVolumeTotal/(iniFluidVolumeTotal+VSMALL)) << " % )" << nl;
    }

    if ( correctVolByHeight > -1)
    {
        label aPatchID = interface.aPatchID();

        scalar volCorr = curFluidVolumePhase1/iniFluidVolumePhase1;

        scalarField volCorrDeltaH = mesh.C().boundaryField()[aPatchID].component(correctVolByHeight)
                                    * (1.0 - volCorr);

        scalarField volCorrPhi = mesh.magSf().boundaryField()[aPatchID]
                                * volCorrDeltaH
                                / runTime.deltaT().value()
                                + fvc::meshPhi(U)().boundaryField()[aPatchID];

        interface.movePoints(volCorrPhi);

        curFluidVolumePhase1 = fvc::domainIntegrate(fluidIndicator).value();
        fluidVolumeCurrentSubDict.set<scalar>("phase1", curFluidVolumePhase1);

        Info << "Fixed volume of phase1!" << endl;

        Info << "Absolute volume of phase1: "
             << curFluidVolumePhase1 << " ( "
             << (100*curFluidVolumePhase1/(iniFluidVolumePhase1+VSMALL)) << " % )" << nl;
    }