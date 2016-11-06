
    tmp<volScalarField> tmuEff(new volScalarField(mu));

    tmuEff() += rho*storage.turbulence().nut();
    tmuEff().rename("muEff");
    tmuEff().correctBoundaryConditions();

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(fvc::interpolate(rho)*phi, U, "div(phi,U)")
      - fvm::laplacian(tmuEff(), U)
      - fvc::div(tmuEff()*dev(fvc::grad(U)().T()))
//       - (fvc::grad(U) & fvc::grad(tmuEff()))
    );

    tmuEff.clear();

    if
    (
        storage.item_F().enabled()
     && settings.UEqnVolumeForce
    )
    {
        if (Control::debug)
        {
            Info << "interTrackApp::Control | UEqn.H : "
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
                - (fvc::snGrad(p)*mesh.magSf())
                )
            );
        }
    }
