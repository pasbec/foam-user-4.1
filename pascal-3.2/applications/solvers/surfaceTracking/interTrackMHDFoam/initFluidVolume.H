
    Info << "Creating/calculating fluid volume(s)" << nl << endl;
    IOdictionary fluidVolumeDict
    (
        IOobject
        (
            "fluidVolume",
            runTime.timeName(),
            "uniform",
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE,
            true
        )
    );

    scalar iniFluidVolumePhase1 = fvc::domainIntegrate(fluidIndicator).value();
    scalar iniFluidVolumePhase2 = fvc::domainIntegrate(1.0-fluidIndicator).value();
    scalar iniFluidVolumeTotal = iniFluidVolumePhase1 + iniFluidVolumePhase2;

    if (fluidVolumeDict.isDict("initial"))
    {
        dictionary fluidVolumeInitialSubDict(fluidVolumeDict.subDict("initial"));

        iniFluidVolumePhase1 = fluidVolumeInitialSubDict.lookupOrAddDefault<scalar>("phase1", iniFluidVolumePhase1);
        iniFluidVolumePhase2 = fluidVolumeInitialSubDict.lookupOrAddDefault<scalar>("phase2", iniFluidVolumePhase2);
        iniFluidVolumeTotal = fluidVolumeInitialSubDict.lookupOrAddDefault<scalar>("total", iniFluidVolumeTotal);
    }
    else
    {
        dictionary fluidVolumeInitialSubDict;

        fluidVolumeInitialSubDict.add<scalar>("phase1", iniFluidVolumePhase1);
        fluidVolumeInitialSubDict.add<scalar>("phase2", iniFluidVolumePhase2);
        fluidVolumeInitialSubDict.add<scalar>("total", iniFluidVolumeTotal);

        fluidVolumeDict.add("initial", fluidVolumeInitialSubDict);
    }


    scalar curFluidVolumePhase1 = iniFluidVolumePhase1;
    scalar curFluidVolumePhase2 = iniFluidVolumePhase2;
    scalar curFluidVolumeTotal = iniFluidVolumeTotal;

    if (fluidVolumeDict.isDict("current"))
    {
        dictionary fluidVolumeCurrentSubDict(fluidVolumeDict.subDict("current"));

        fluidVolumeCurrentSubDict.add<scalar>("phase1", curFluidVolumePhase1, true);
        fluidVolumeCurrentSubDict.add<scalar>("phase2", curFluidVolumePhase2, true);
        fluidVolumeCurrentSubDict.add<scalar>("total", curFluidVolumeTotal, true);
    }
    else
    {
        dictionary fluidVolumeCurrentSubDict;

        fluidVolumeCurrentSubDict.add<scalar>("phase1", curFluidVolumePhase1);
        fluidVolumeCurrentSubDict.add<scalar>("phase2", curFluidVolumePhase2);
        fluidVolumeCurrentSubDict.add<scalar>("total", curFluidVolumeTotal);

        fluidVolumeDict.add("current", fluidVolumeCurrentSubDict);
    }

    dictionary& fluidVolumeCurrentSubDict = fluidVolumeDict.subDict("current");