
    tmp<volScalarField> tmu = mu;

    if(storage.item_turbulence().enabled())
    {
        tmu() += rho*storage.turbulence().nut();
        tmu().rename("muEff");
    }

    tmp<fvVectorMatrix> UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(fvc::interpolate(rho)*phi, U, "div(phi,U)")
      - fvm::laplacian(tmu(), U)
      - fvc::div(tmu()*dev(fvc::grad(U)().T()))
//       - (fvc::grad(U) & fvc::grad(tmu()))
    );

    tmu.clear();

    if
    (
        storage.item_F().enabled()
     && settings.UEqnVolumeForce
    )
    {
        if (Manager::debug)
        {
            Info << "interTrackManager | UEqn.H : "
                << "Adding source term for volume force F in UEqn."
                << endl;
        }

        UEqn() -= storage.F();
    }

    UEqn().relax();

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
            solve(UEqn() == -fvc::grad(p));
        }
    }

    if (settings.UpCoupled)
    {
        // Init block matrix
        UpEqn.set(new fvVector4Matrix(storage.Up()));

        // Add UEqn into block matrix
        UpEqn().insertEquation(0, UEqn());
    }
