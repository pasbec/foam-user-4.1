
    // fvSolution
    dictionary& fvSolution = mesh.solutionDict();

    // AV dictionary
    dictionary AVdict = fvSolution.subDict("AV");

    // AV method
    word AVmethod = AVdict.lookupOrDefault<word>("method", "coupled");

    dictionary AVmethodDict = AVdict.subDict(AVmethod);
