
    // AV dictionary
    dictionary& AVdict = mesh.solutionDict().subDict("AV");

    // AV settings
    int minIter = AVdict.lookupOrDefault<int>("minIter", 0);
    int maxIter = AVdict.lookupOrDefault<int>("maxIter", 100);
    scalar tol = AVdict.lookupOrDefault<scalar>("tolerance", 1e-04);

    // Copy local solver dictionaries
    dictionary Adict = mesh.solutionDict().subDict("solvers").subDict("A");
    dictionary Vdict = cMesh.solutionDict().subDict("solvers").subDict("V");

    // Init current tolerances
    scalar curTol = 1.0;
    scalar curTolA = curTol - SMALL;
    scalar curTolV = curTol - SMALL;

    // AV iterations
    for (label iter = 0; iter < maxIter; iter++)
    {
        if ((iter >= minIter) && (curTol <= tol)) break;

        // TODO; curTolA
        if (curTolV < curTol)
        {
            curTol /= 10.0;

            Info << nl << "Tolerance = " << curTol
                << " / " << tol << endl << nl;

            Adict.set<scalar>("tolerance", curTol);
            Vdict.set<scalar>("tolerance", curTol);
        }

        // Solve for A
#       include "AEqn.H"

        // Solve for V
#       include "VEqn.H"
    }
