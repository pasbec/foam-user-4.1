
    // fvSolution
    dictionary& fvSolution = mesh.solutionDict();

    // AV dictionary
    dictionary AVdict = fvSolution.subDict("AV");

    // AV method
    word AVmethod = AVdict.lookupOrDefault<word>("method", "coupled");

    dictionary AVmethodDict = AVdict.subDict(AVmethod);

    // Read scale factors
    scalar scaleA = AVmethodDict.lookupOrDefault<scalar>("scaleA", 1.0);
    scalar scaleV = AVmethodDict.lookupOrDefault<scalar>("scaleV", 1.0);
