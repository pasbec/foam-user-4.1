
    // Momentum equation
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U, "ddt(U)")
      + rho*fvm::div(phi, U, "div(phi,U)")
      + rho*turbulence->divDevReff(U)
      - rho*FL * lorentzForceVolumeFactor
      - fvc::grad(pB) * lorentzForceVolumeFactor * lorentzForceRotationalFactor
    );

    if (outerCorr != nOuterCorr - 1)
    {
        Info << "Relaxing UEqn with factor: " << mesh.solutionDict().relaxationFactor("U") << nl;
        UEqn.relax(mesh.solutionDict().relaxationFactor("U"));
    }

    blockVectorMatrix pInU(fvm::grad(p,"BLOCK[fvm::grad(p)]"));