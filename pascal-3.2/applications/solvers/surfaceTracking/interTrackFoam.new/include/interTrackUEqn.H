
    tmp<fvVectorMatrix> UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(fvc::interpolate(rho)*phi, U, "div(phi,U)")
//       - fvm::laplacian(mu, U)
// TODO: Turbulence
      + rho*turbulence.divDevReff(U)
// TODO: Check the way turbulence is introduced in the U-equation
//       - fvm::laplacian(nuEff(), U)
//       - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );

    UEqn().relax();

    if (control.momentumPredictor())
    {
        solve(UEqn() == -fvc::grad(p));
    }
