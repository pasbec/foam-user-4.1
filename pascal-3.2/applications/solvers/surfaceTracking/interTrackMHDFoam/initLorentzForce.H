
    Info << "Creating/calculating Lorentz force control variables" << nl << endl;
    IOdictionary lorentzForceTimeDict
    (
        IOobject
        (
            "lorentzForce",
            runTime.timeName(),
            "uniform",
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE,
            true
        )
    );

    // Init output index for Lorentz force update
    int lorentzForceUpdateOutputTimeIndex = lorentzForceTimeDict.lookupOrAddDefault<int>("outputTimeIndex", 0);

    // Init last time of lorentz force update time cycle
    scalar lorentzForceUpdateTimeCycleLastUpdateTime = lorentzForceTimeDict.lookupOrAddDefault<scalar>("lastUpdateTime", runTime.value());