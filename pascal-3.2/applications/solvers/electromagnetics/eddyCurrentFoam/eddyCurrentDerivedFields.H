
    // Derived fields in base region
    {
        setRegionScope(control.base());

        // Magnetic field density
        BRe == fvc::curl(ARe);
        BIm == fvc::curl(AIm);

        // Eddy current density
        jRe ==   control.omega() * sigma * AIm - sigma * VReGrad;
        jIm == - control.omega() * sigma * ARe - sigma * VImGrad;

        // Time-averaged Lorentz-force
        FL == 0.5 * ( (jRe ^ BRe) + (jIm ^ BIm) );

        // Time-averaged magnetic pressure
        pB == 0.5 * physicalConstant::rMu0
            * 0.5 * ( (BRe & BRe) + (BIm & BIm) );
    }

    // Derived fields in conductor region
    {
        FL_.mapExtrapolate(control.conductor());
        pB_.mapExtrapolate(control.conductor());
    }
