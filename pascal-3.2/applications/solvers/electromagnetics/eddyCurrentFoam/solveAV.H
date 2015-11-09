
    // Copy local solver dictionaries
    dictionary Adict = mesh.solutionDict().subDict("solvers").subDict("A");
    dictionary Vdict = cMesh.solutionDict().subDict("solvers").subDict("V");

    // Reset current tolerance
    curTol = 1.0;

    // AV iterations
    for (int iter = 0; iter < maxIter; iter++)
    {
        if ((iter >= minIter) && (res < tol)) break;

        if (res < curTol)
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

        // Update residual
        res = max(resA, resV);
    }
