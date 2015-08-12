
    IOdictionary mfTimeDict
    (
        IOobject
        (
            "magnetic",
            runTime.timeName(),
            "uniform",
            mesh_[defaultRegionID],
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE,
            true
        )
    );

    // Init runtime update switch
    bool mfUpdate = true;

    // Init update counter
    int mfUpdateCounter =
        mfTimeDict.lookupOrAddDefault<int>("updateCounter", 0);

    // Init output index for magnetic update
    // We need this as there is currently no access function to
    // private variable outputTimeIndex_ of TimeState
    int mfUpdateOutputTimeIndex =
        mfTimeDict.lookupOrAddDefault<int>("outputTimeIndex", 0);

    // Init last time of magnetic update time cycle
    scalar mfUpdateLastUpdateTime =
        mfTimeDict.lookupOrAddDefault<scalar>("lastUpdateTime", -VGREAT);
