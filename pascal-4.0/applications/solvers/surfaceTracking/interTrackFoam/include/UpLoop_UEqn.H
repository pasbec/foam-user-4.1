
    surfaceScalarField muEff
    (
        "muEff",
        transport.muf()
      + fvc::interpolate(rho*turbulence.nut())
    );

    surfaceScalarField rhoPhi
    (
        "rho*phi",
        fvc::interpolate(rho)*phi
    );

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      - fvm::laplacian(muEff, U)
      - (fvc::grad(U) & fvc::grad(muEff))
    //- fvc::div(muEff*(fvc::interpolate(dev(fvc::grad(U))) & mesh.Sf()))
    );

    UEqn.relax();

    if (control.momentumPredictor())
    {
        if (storage.item_F().enabled())
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
