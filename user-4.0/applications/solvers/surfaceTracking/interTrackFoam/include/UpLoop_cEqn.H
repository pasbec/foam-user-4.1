
    if (storage.item_c().enabled())
    {
        uniformDimensionedScalarField& Dc = storage.Dc();
        volScalarField& c = storage.c();

        fvScalarMatrix cEqn
        (
            fvm::ddt(c)
          + fvm::div(phi, c)
          - fvm::laplacian(Dc, c)
        );

        cEqn.relax();

        cEqn.solve();

        if (settings.relToUinf)
        {
            // WARNING: DO NOT TOUCH boundary values!
            c.internalField() -= max(c).value();
            c.correctBoundaryConditions();
        }
    }
