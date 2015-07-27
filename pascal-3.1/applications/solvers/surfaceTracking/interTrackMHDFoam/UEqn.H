
    rhof = fvc::interpolate(rho);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U, "ddt(U)")
      + fvm::div(rhof*phi, U, "div(phi,U)")
      + rho*turbulence->divDevReff(U)
    );

    if (lorentzForceNonAlberto)
    {
        UEqn = UEqn
             - rho*FL * lorentzForceVolumeFactor
             - fvc::grad(pB) * lorentzForceVolumeFactor * lorentzForceRotationalFactor;
    }

    word USolverDict = "UFinal";

    if (outerCorr != nOuterCorr - 1)
    {
        USolverDict = "U";

        Info << "Relaxing UEqn with factor: " << mesh.solutionDict().relaxationFactor("U") << nl;
        UEqn.relax(mesh.solutionDict().relaxationFactor("U"));
    }

    if (lorentzForceNonAlberto)
    {
        solve
        (
            UEqn == -fvc::grad(p),
            mesh.solutionDict().solver(USolverDict)
        );
    }
    else
    {
        solve
        (
            UEqn == fvc::reconstruct
            (
                - ( fvc::snGrad(p) * mesh.magSf() )
                + ( fvc::interpolate(rho*FL) & mesh.Sf() ) * lorentzForceVolumeFactor
                + ( fvc::snGrad(pB) * mesh.magSf() ) * lorentzForceVolumeFactor * lorentzForceRotationalFactor
            ),
            mesh.solutionDict().solver(USolverDict)
        );
    }

    UcorrByrAU.storePrevIter();