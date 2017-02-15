
    // Solve A0 using Biot-Savart law
    if (globalSettings.biotSavart)
    {
        SM_REGIONSCOPE(CONDUCTOR);

        if (Control::debug)
        {
            Info<< Control::typeName << " | A0BiotSavart.H : "
                << "Calculate impressed vector potential"
                << " via Biot-Savart's law."
                << endl;
        }

        // Construct Biot-Savart class
        edgeBiotSavart eBiotSavart(globalMesh[CONDUCTOR]);
//         eBiotSavart.writeVTK();

        // Calculate A0 at boundaries
        globalStorage.A0Re()[CONDUCTOR] == eBiotSavart.A(edgeBiotSavart::REAL, true);
        globalStorage.A0Im()[CONDUCTOR] == eBiotSavart.A(edgeBiotSavart::IMAG, true);
    }
