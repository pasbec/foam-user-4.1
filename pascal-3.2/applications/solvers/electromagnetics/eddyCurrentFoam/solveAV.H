
    // AV iterations
    int minIter = AVdict.lookupOrDefault<int>("minIter", 0);
    int maxIter = AVdict.lookupOrDefault<int>("maxIter", 100);
    scalar tolerance = AVdict.lookupOrDefault<scalar>("tolerance", 1e-04);

    // Tolerance

    for (label iter = 0; iter < maxIter; iter++)
    {
        // Solve for A
#       include "AEqn.H"

        // Solve for V
#       include "VEqn.H"
    }
