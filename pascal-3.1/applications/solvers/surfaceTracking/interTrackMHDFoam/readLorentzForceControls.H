
    dictionary& lorentzForceSolutionDict(pimple.subDict("lorentzForce"));

    scalar lorentzForceRotationalFactor =
        lorentzForceSolutionDict.lookupOrDefault<scalar>("rotationalFactor", 0.0);

    scalar lorentzForceVolumeFactor =
        lorentzForceSolutionDict.lookupOrDefault<scalar>("volumeFactor", 1.0);

    Switch lorentzForceNonAlberto =
        lorentzForceSolutionDict.lookupOrDefault<Switch>("nonAlberto", true);


    Switch lorentzForceUpdateSwitch =
        lorentzForceSolutionDict.lookupOrDefault<Switch>("updateSwitch", false);

    // HINT: Maximum unsigned 32bit int number is: 2147483647
    int lorentzForceUpdateOutputCycle =
        lorentzForceSolutionDict.lookupOrDefault<int>("updateOutputCycle", std::numeric_limits<int>::max());

    int lorentzForceUpdateTimeStepCycle =
        lorentzForceSolutionDict.lookupOrDefault<int>("updateTimeStepCycle", std::numeric_limits<int>::max());

    scalar lorentzForceUpdateTimeCycle =
        lorentzForceSolutionDict.lookupOrDefault<scalar>("updateTimeCycle", VGREAT);

    string lorentzForceUpdateUtility =
        lorentzForceSolutionDict.lookup("updateUtility");