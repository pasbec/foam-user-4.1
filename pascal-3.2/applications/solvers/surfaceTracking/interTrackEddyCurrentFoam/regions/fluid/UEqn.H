
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U, "ddt(U)")
      + rho*fvm::div(phi, U, "div(phi,U)")
      + rho*turbulence->divDevReff(U)
    );

    if (methodUp == "default")
    {
        UEqn -= F;
    }

    word USolverDict = "UFinal";

    if (outerCorr != nOuterCorr - 1)
    {
        USolverDict = "U";

        scalar alpha = mesh.solutionDict().relaxationFactor("U");

        if (alpha < 1.0)
        {
            Info << "Relaxing U with factor: "
                << alpha << nl;

            UEqn.relax(alpha);
        }
    }

    if (methodUp == "default")
    {
        solve
        (
            UEqn == -fvc::grad(p),
            mesh.solutionDict().solver(USolverDict)
        );
    }
    else if (methodUp == "alberto")
    {
        solve
        (
            UEqn == fvc::reconstruct
            (
                - ( fvc::snGrad(p) * mesh.magSf() )
                + ( fvc::interpolate(F) & mesh.Sf() )
            ),
            mesh.solutionDict().solver(USolverDict)
        );
    }
