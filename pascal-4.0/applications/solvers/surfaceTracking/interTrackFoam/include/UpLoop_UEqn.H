
    volScalarField mu
    (
        "mu",
        rho*transport.nu()
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

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      - fvm::laplacian(muEff, U)
      - (fvc::grad(U) & fvc::grad(muEff))
    //- fvc::div(muEff*(fvc::interpolate(dev(fvc::grad(U))) & mesh.Sf()))
    );

    if (!settings.diffusiveFlow)
    {

// TODO: Two fluids: Is this ok?
        surfaceScalarField rhoPhi
        (
            "rho*phi",
            fvc::interpolate(rho)*phi
        );

        UEqn += fvm::div(rhoPhi, U);
    }

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
