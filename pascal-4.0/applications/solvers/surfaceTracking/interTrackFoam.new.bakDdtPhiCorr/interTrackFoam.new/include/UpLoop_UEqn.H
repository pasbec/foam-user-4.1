
    // Time derivative matrix
    fvVectorMatrix ddtUEqn(fvm::ddt(rho, U));

    tmp<volScalarField> tmuEff(new volScalarField(mu));

    if(storage.item_turbulence().enabled())
    {
        tmuEff() += rho*storage.turbulence().nut();
        tmuEff().rename("muEff");
    }

    tmuEff().correctBoundaryConditions();

    // Convection-diffusion matrix
    tmp<fvVectorMatrix> HUEqn
    (
        fvm::div(fvc::interpolate(rho)*phiHbyA, U, "div(phi,U)")
      - fvm::laplacian(tmuEff(), U)
      - fvc::div(tmuEff()*dev(fvc::grad(U)().T()))
//       - (fvc::grad(U) & fvc::grad(tmuEff()))
    );

    tmuEff.clear();

//     if
//     (
//         storage.item_F().enabled()
//      && settings.UEqnVolumeForce
//     )
//     {
//         if (Control::debug)
//         {
//             Info << "interTrackApp::Control | UEqn.H : "
//                 << "Adding source term for volume force F in UEqn."
//                 << endl;
//         }
//
//         HUEqn() -= storage.F();
//     }

// TODO: Relaxation
//     UEqn().relax();

// TODO: Predictor step
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
//                 ddtUEqn + HUEqn() == fvc::reconstruct
//                 (
//                 - (fvc::snGrad(p)*mesh.magSf())
//                 + (fvc::interpolate(storage.F()) & mesh.Sf())
//                 )
//             );
//         }
//         else
//         {
//             solve(ddtUEqn + HUEqn() == -fvc::grad(p));
//         }
//     }

// TODO: Coupled solution
//     if (settings.UpCoupled)
//     {
//         // Init block matrix
//         UpEqn.set(new fvVector4Matrix(storage.Up()));
//
//         // Add UEqn into block matrix
//         UpEqn().insertEquation(0, UEqn());
//     }
