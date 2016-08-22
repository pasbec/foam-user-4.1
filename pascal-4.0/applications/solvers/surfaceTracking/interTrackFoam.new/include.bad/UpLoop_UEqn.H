
    tmp<volScalarField> tmuEff(new volScalarField(mu));

    if(storage.item_turbulence().enabled())
    {
        tmuEff() += rho*storage.turbulence().nut();
        tmuEff().rename("muEff");
    }

    tmuEff().correctBoundaryConditions();

    // Time derivative matrix of U-equation
    tmp<fvVectorMatrix> ddtUEqn(fvm::ddt(rho, U));

    // Matrix of U-equation without time derivative
    tmp<fvVectorMatrix> HUEqn
    (
        fvm::div(fvc::interpolate(rho)*phi, U, "div(phi,U)")
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
        if (Manager::debug)
        {
            Info << "interTrackManager | UEqn.H : "
                << "Adding source term for volume force F in HUEqn."
                << endl;
        }

        HUEqn() -= storage.F();
    }

// TODO: REMOVE ddtPhiCorr
    // Under-relaxation factor of U-equation
//     scalar UUrf =
//         mesh.solutionDict().equationRelaxationFactor(U.select(control.finalIter()));
    scalar UUrf = 1.0;

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
                ddtUEqn()
              + relax(HUEqn(), UUrf)
             == fvc::reconstruct
                (
                - (fvc::snGrad(p)*mesh.magSf())
                + (fvc::interpolate(storage.F()) & mesh.Sf())
                ),
                mesh.solutionDict().solver((U.select(control.finalIter())))
            );
        }
        else
        {
            solve
            (
                ddtUEqn()
              + relax(HUEqn(), UUrf)
             == fvc::reconstruct
                (
                - (fvc::snGrad(p)*mesh.magSf())
                ),
                mesh.solutionDict().solver(U.select(control.finalIter()))
            );
        }
    }

// TODO: REMOVE ddtPhiCorr
//     if (settings.UpCoupled)
//     {
//         // Init block matrix
//         UpEqn.set(new fvVector4Matrix(storage.Up()));
//
//         // Add UEqn into block matrix
//         UpEqn().insertEquation(0, UEqn());
//     }
