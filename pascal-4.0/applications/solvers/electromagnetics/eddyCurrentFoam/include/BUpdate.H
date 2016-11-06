
    // Calculate B
    {
        if (Control::debug)
        {
            Info<< Control::typeName << " | BUpdate.H : "
                << "Calculate magnetic field."
                << endl;
        }

        globalStorage.BRe()[DEFAULT] = fvc::curl(globalStorage.ARe()[DEFAULT]);
        globalStorage.BIm()[DEFAULT] = fvc::curl(globalStorage.AIm()[DEFAULT]);

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

            // Update boundary field of B
            globalStorage.BRe().extrapolateBoundaryField(CONDUCTOR);
            globalStorage.BIm().extrapolateBoundaryField(CONDUCTOR);
        }
    }
