
    if (storage.item_c().enabled())
    {
        volScalarField& c = storage.c();

        const dictionary& concentrationDict =
            transport.subDict("concentration");

        dimensionedScalar betaRho(concentrationDict.lookup("betaRho"));
        dimensionedScalar cRefRho(concentrationDict.lookup("cRefRho"));

// TODO: Two fluids: Use fluidIndicator!
        rho = transport.rho() + betaRho*(c - cRefRho);
    }
