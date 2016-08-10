
    // Solve A0 using Biot-Savart law
    if (globalSettings.biotSavart)
    {
        // Construct Biot-Savart class
        edgeBiotSavart eBiotSavart(globalMesh[CONDUCTOR]);

        eBiotSavart.writeVTK();

        volVectorField& A0Re = globalStorage.A0Re()[CONDUCTOR];
        volVectorField& A0Im = globalStorage.A0Im()[CONDUCTOR];

        // Calculate A0
        A0Re = eBiotSavart.A(A0Re.name(), edgeBiotSavart::REAL);
        A0Im = eBiotSavart.A(A0Im.name(), edgeBiotSavart::IMAG);

        // Reverse map A0 to default region
        globalStorage.A0Re().rmap(CONDUCTOR);
        globalStorage.A0Im().rmap(CONDUCTOR);
    }
