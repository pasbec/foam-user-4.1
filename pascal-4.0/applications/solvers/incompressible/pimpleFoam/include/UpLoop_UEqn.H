    // Solve the momentum equation

    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + turbulence.divDevReff()
    );

    UEqn.relax();

    if (control.momentumPredictor())
    {
        solve
        (
            UEqn == fvc::reconstruct
            (
              - (fvc::snGrad(p)*mesh.magSf())
              + (fvc::interpolate(storage.F())/storage.rho() & mesh.Sf())
            )
        );
    }
