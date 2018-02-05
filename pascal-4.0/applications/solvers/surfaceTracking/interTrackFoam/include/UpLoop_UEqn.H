
    volScalarField mu
    (
        "mu",
        rho*turbulence.nu()
    );

    surfaceScalarField muf
    (
        "muf",
        fvc::interpolate(mu)
    );

    surfaceScalarField muEff
    (
        "muEff",
        muf
      + fvc::interpolate(rho*turbulence.nut())
    );

// TODO: Two fluids: Is this ok?
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

    // Hele-Shaw Poisson drag force
    if (storage.item_heleShawGapWidth().enabled())
    {
        if (Control::debug)
        {
            Info<< Control::typeName << " | UEqn.H : "
                << "Adding Hele-Shaw Poisson drag term to UEqn."
                << endl;
        }

        uniformDimensionedScalarField& HSGW = storage.heleShawGapWidth();

        UEqn += fvm::SuSp(3.0/pow(HSGW, 2)*mu, U);
    }

    UEqn.relax();

    if (control.momentumPredictor())
    {
        if (storage.item_F().enabled())
        {
            solve
            (
                UEqn == fvc::reconstruct
                (
                    (
                      - gxf*fvc::snGrad(rho)
                      - fvc::snGrad(p)
                    )*mesh.magSf()
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
                    (
                      - gxf*fvc::snGrad(rho)
                      - fvc::snGrad(p)
                    )*mesh.magSf()
                )
            );
        }
    }
