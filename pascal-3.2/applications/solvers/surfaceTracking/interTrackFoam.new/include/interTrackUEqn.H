
    tmp<volScalarField> tmu = mu;

    if(storage.is_turbulence())
    {
        tmu() += rho * storage.turbulence().nut();
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

    if(storage.is_F())
    {
        UEqn() -= storage.F();
    }

    UEqn().relax();

    if (control.momentumPredictor())
    {
        solve(UEqn() == -fvc::grad(p));
    }
