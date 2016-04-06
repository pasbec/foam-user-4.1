
    // Derived fields in base region
    {
        // Magnetic field density
        BRe == fvc::curl(ARe);
        BIm == fvc::curl(AIm);

        // Eddy current density
        jRe ==   omega0 * sigma * AIm - sigma * VReGrad;
        jIm == - omega0 * sigma * ARe - sigma * VImGrad;

        // Time-averaged Lorentz-force
        FL == 0.5 * ( (jRe ^ BRe) + (jIm ^ BIm) );

        // Time-averaged magnetic pressure
        pB == 0.5 * physicalConstant::rMu0
            * 0.5 * ( (BRe & BRe) + (BIm & BIm) );
    }

    // Derived fields in conductor region
    {
        FL.mapExtrapolate(control.conductor());
        pB.mapExtrapolate(control.conductor());
// TODO TEST FIXME Extrapolation sometimes fails
//         FL.mapCopyInternal(control.conductor());
//         pB.mapCopyInternal(control.conductor());
    }
