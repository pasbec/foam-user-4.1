
    tmp<volScalarField> muEff = rho * turbulence.nuEff();
    muEff().rename("muEff");

    tmp<fvVectorMatrix> UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(fvc::interpolate(rho)*phi, U, "div(phi,U)")
      - fvm::laplacian(muEff(), U)
      - fvc::div(muEff()*dev(fvc::grad(U)().T()))
//       - (fvc::grad(U) & fvc::grad(muEff()))
      - F
    );

    muEff.clear();

    UEqn().relax();

    if (control.momentumPredictor())
    {
        solve(UEqn() == -fvc::grad(p));
    }
