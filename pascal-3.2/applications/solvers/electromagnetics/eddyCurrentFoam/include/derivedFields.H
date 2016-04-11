
    // Derived fields in base region
    {
        SM_REGIONSCOPE(DEFAULT);

        volScalarField& sigma = globalStorage.sigma()[DEFAULT];

        volVectorField& ARe = globalStorage.ARe()[DEFAULT];
        volVectorField& AIm = globalStorage.AIm()[DEFAULT];

        volVectorField& VReGrad = globalStorage.VReGrad()[DEFAULT];
        volVectorField& VImGrad = globalStorage.VImGrad()[DEFAULT];

        volVectorField& FL = globalStorage.FL()[DEFAULT];
        volScalarField& pB = globalStorage.pB()[DEFAULT];

        // Magnetic field density
        storage.BRe() == fvc::curl(ARe);
        storage.BIm() == fvc::curl(AIm);

        // Eddy current density
        storage.jRe() ==   omega0 * sigma * AIm - sigma * VReGrad;
        storage.jIm() == - omega0 * sigma * ARe - sigma * VImGrad;

        // Time-averaged Lorentz-force
        FL == 0.5 * ( (storage.jRe() ^ storage.BRe())
                    + (storage.jIm() ^ storage.BIm()) );

        // Time-averaged magnetic pressure
        pB == 0.5 * physicalConstant::rMu0
            * 0.5 * ( (storage.BRe() & storage.BRe())
                    + (storage.BIm() & storage.BIm()) );
    }

    // Derived fields in conductor region
    {
        globalStorage.FL().mapExtrapolate(CONDUCTOR);
        globalStorage.pB().mapExtrapolate(CONDUCTOR);
// // TODO FIXME Extrapolation sometimes fails
//         globalStorage.FL().mapCopyInternal(CONDUCTOR);
//         globalStorage.pB().mapCopyInternal(CONDUCTOR);
    }
