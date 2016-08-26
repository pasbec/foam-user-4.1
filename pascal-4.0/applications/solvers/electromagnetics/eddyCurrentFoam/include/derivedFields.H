
    // Derived fields in conductor region
    {
        SM_REGIONSCOPE(CONDUCTOR);

        volScalarField& sigma = globalStorage.sigma()[CONDUCTOR];

        volVectorField& ARe = globalStorage.ARe()[CONDUCTOR];
        volVectorField& AIm = globalStorage.AIm()[CONDUCTOR];

        volVectorField& FL = globalStorage.FL()[CONDUCTOR];
        volScalarField& pB = globalStorage.pB()[CONDUCTOR];

        volVectorField& BRe = globalStorage.BRe()[CONDUCTOR];
        volVectorField& BIm = globalStorage.BIm()[CONDUCTOR];

        // Eddy current density
        storage.jRe() =   omega0 * sigma * AIm;
        storage.jIm() = - omega0 * sigma * ARe;

        if (control.meshIs3D())
        {
            storage.jRe() -= sigma * globalStorage.VReGrad()[CONDUCTOR];
            storage.jIm() -= sigma * globalStorage.VImGrad()[CONDUCTOR];
        }

        if (globalSettings.biotSavart)
        {
            storage.jRe() += omega0 * sigma * globalStorage.A0Im()[CONDUCTOR];
            storage.jIm() -= omega0 * sigma * globalStorage.A0Re()[CONDUCTOR];
        }

        // Time-averaged Lorentz-force
        FL = 0.5 * ( (storage.jRe() ^ BRe)
                   + (storage.jIm() ^ BIm) );

        // Time-averaged magnetic pressure
        pB = 0.5 * physicalConstant::rMu0
           * 0.5 * ( (BRe & BRe)
                   + (BIm & BIm) );
    }

    // Derived fields in default region
    {
        globalStorage.FL().rmap(CONDUCTOR);
        globalStorage.pB().rmap(CONDUCTOR);
    }

