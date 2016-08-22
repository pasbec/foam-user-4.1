
    tmp<volScalarField> tmuEff(new volScalarField(mu));

    if(storage.item_turbulence().enabled())
    {
        tmuEff() += rho*storage.turbulence().nut();
        tmuEff().rename("muEff");
    }

    tmuEff().correctBoundaryConditions();

    // Time derivative matrix of U-equation
    fvVectorMatrix ddtUEqn(fvm::ddt(rho, U));

    // Matrix of U-equation without time derivative
    fvVectorMatrix HUEqn
    (
        fvm::div(fvc::interpolate(rho)*phiH, U, "div(phi,U)")
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

        HUEqn -= storage.F();
    }

// TODO: REMOVE ddtPhiCorr
//     HUEqn.relax();
//
//     if (control.momentumPredictor())
//     {
//         if
//         (
//             storage.item_F().enabled()
//         && !settings.UEqnVolumeForce
//         )
//         {
//             solve
//             (
//                 ddtUEqn + HUEqn == fvc::reconstruct
//                 (
//                 - (fvc::snGrad(p)*mesh.magSf())
//                 + (fvc::interpolate(storage.F()) & mesh.Sf())
//                 )
//             );
//         }
//         else
//         {
//             solve(ddtUEqn + HUEqn == -fvc::grad(p));
//         }
//     }

// TODO: REMOVE ddtPhiCorr
//     if (settings.UpCoupled)
//     {
//         // Init block matrix
//         UpEqn.set(new fvVector4Matrix(storage.Up()));
//
//         // Add UEqn into block matrix
//         UpEqn().insertEquation(0, UEqn());
//     }
