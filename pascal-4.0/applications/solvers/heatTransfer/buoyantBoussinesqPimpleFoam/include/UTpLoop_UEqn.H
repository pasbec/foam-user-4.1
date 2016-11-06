    // Solve the momentum equation

    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + turbulence.divDevReff()
    );

    if
    (
        storage.item_F().enabled()
     && settings.UEqnVolumeForce
    )
    {
        if (Control::debug)
        {
            Info<< Control::typeName << " | UEqn.H : "
                << "Adding source term for volume force F in UEqn."
                << endl;
        }

        UEqn -= storage.F();
    }

    UEqn.relax();

    if (control.momentumPredictor())
    {
        if
        (
            storage.item_F().enabled()
        && !settings.UEqnVolumeForce
        )
        {
            solve
            (
                UEqn == fvc::reconstruct
                (
                  (fvc::interpolate(rhok)*(g & mesh.Sf()))
                - (fvc::snGrad(p)*mesh.magSf())
                + (fvc::interpolate(storage.F()) & mesh.Sf())
                )
            );
        }
        else
        {
            solve
            (
                UEqn == fvc::reconstruct
                (
                  (fvc::interpolate(rhok)*(g & mesh.Sf()))
                - (fvc::snGrad(p)*mesh.magSf())
                )
            );
        }
    }
