
    if (storage.item_c().enabled())
    {
        uniformDimensionedScalarField& Dc = storage.Dc();
        volScalarField& c = storage.c();

        fvScalarMatrix cEqn
        (
            fvm::ddt(c)             // Unsteady term
          + fvm::div(phi, c)        // Convective term
          - fvm::laplacian(Dc, c)   // Diffusive term
        );

        cEqn.relax();

        cEqn.solve();
    }
