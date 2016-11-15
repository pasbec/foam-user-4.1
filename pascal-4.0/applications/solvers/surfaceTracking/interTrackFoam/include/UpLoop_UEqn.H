
    tmp<volScalarField> tmuEff(new volScalarField(mu));

    tmuEff() += rho*turbulence.nut();
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

    UEqn.relax();

    if (control.momentumPredictor())
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
