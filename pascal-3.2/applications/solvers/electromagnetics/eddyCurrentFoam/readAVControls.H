
    // AV dictionary
    dictionary& AVdict = mesh.solutionDict().subDict("AV");

    // AV settings
    int minIter = AVdict.lookupOrDefault<int>("minIter", 0);
    int maxIter = AVdict.lookupOrDefault<int>("maxIter", 100);
    scalar tol = AVdict.lookupOrDefault<scalar>("tolerance", 1e-04);
