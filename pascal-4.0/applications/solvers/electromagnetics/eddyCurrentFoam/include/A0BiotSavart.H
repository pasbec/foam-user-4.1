
    // Solve A0 using Biot-Savart law
    if (globalSettings.biotSavart)
    {
        if (Control::debug)
        {
            Info << "eddyCurrentApp::Control | A0BiotSavart.H : "
                << "Calculate impressed vector potential"
                << " via Biot-Savart's law."
                << endl;
        }

// TODO: Put this into manager! and write only once
        // Construct Biot-Savart class
        edgeBiotSavart eBiotSavart(globalMesh[CONDUCTOR]);
//         eBiotSavart.writeVTK();

        // Calculate A0
        eBiotSavart.A(globalStorage.A0Re()[CONDUCTOR], edgeBiotSavart::REAL);
        eBiotSavart.A(globalStorage.A0Im()[CONDUCTOR], edgeBiotSavart::IMAG);

        // Reverse map A0 to default region
        globalStorage.A0Re().rmap(CONDUCTOR);
        globalStorage.A0Im().rmap(CONDUCTOR);
    }
