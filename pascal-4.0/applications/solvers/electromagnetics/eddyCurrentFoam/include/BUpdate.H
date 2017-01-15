
    // Calculate B
    {
        if (Control::debug)
        {
            Info<< Control::typeName << " | BUpdate.H : "
                << "Calculate B."
                << endl;
        }

        globalStorage.BRe()[DEFAULT] = fvc::curl(globalStorage.ARe()[DEFAULT]);
        globalStorage.BIm()[DEFAULT] = fvc::curl(globalStorage.AIm()[DEFAULT]);

// volScalarField divBRe = fvc::div(globalStorage.BRe()[DEFAULT]);
// volScalarField divBIm = fvc::div(globalStorage.BIm()[DEFAULT]);
// divBRe.write();
// divBIm.write();

        // Update B in conductor region
        globalStorage.BRe().map(CONDUCTOR);
        globalStorage.BIm().map(CONDUCTOR);

        // If Biot-Savart IS used, A0 in the non-cunducting region
        // is NOT calculated. Thus, we are not allowed to do fvc::curl(A + A0)
        // in the whole default region.
        if (globalSettings.biotSavart)
        {
            globalStorage.BRe()[CONDUCTOR] += fvc::curl(globalStorage.A0Re()[CONDUCTOR]);
            globalStorage.BIm()[CONDUCTOR] += fvc::curl(globalStorage.A0Im()[CONDUCTOR]);
        }

        // Update boundary field of B
        globalStorage.BRe().extrapolateBoundaryField(CONDUCTOR);
        globalStorage.BIm().extrapolateBoundaryField(CONDUCTOR);

// volScalarField divBReC = fvc::div(globalStorage.BRe()[CONDUCTOR]);
// volScalarField divBImC = fvc::div(globalStorage.BIm()[CONDUCTOR]);
// divBReC.write();
// divBImC.write();
    }
