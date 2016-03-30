
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(fvc::interpolate(rho)*phi, U, "div(phi,U)")
      - fvm::laplacian(mu, U)
// TODO: Turbulence
// TODO: Check the way turbulence is introduced in the U-equation
//       + rho*turbulence->divDevReff(U)
// OR
//       - fvm::laplacian(nuEff(), U)
//       - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );

    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));
    }
