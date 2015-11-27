
    dictionary& magnetic = mesh_[control.base()].solutionDict().subDict("Magnetic");

    // Update settings
    dictionary& mfUpdateDict(magnetic.subDict("update"));

        bool mfUpdateSwitch =
            mfUpdateDict.lookupOrDefault<bool>("active", true);

        int mfUpdateOutputCycle =
            mfUpdateDict.lookupOrDefault<int>("outputCycle", 2147483647);

        int mfUpdateTimeStepCycle =
            mfUpdateDict.lookupOrDefault<int>("timeStepCycle", 2147483647);

        scalar mfUpdateTimeCycle =
            mfUpdateDict.lookupOrDefault<scalar>("timeCycle", VGREAT);
