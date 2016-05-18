
    // Solve A0 using Biot-Savart law
    if (globalSettings.biotSavart)
    {
        // Construct Biot-Savart class
        edgeBiotSavart eBiotSavart(globalMesh[CONDUCTOR]);

        volVectorField& A0Re = globalStorage.A0Re()[CONDUCTOR];
        volVectorField& A0Im = globalStorage.A0Im()[CONDUCTOR];

        // Calculate A0 on boundary (interface) patches
        A0Re.boundaryField()[control.interfaceLabel()] ==
            eBiotSavart.boundaryPatchA
            (
                control.interfaceLabel(),
                edgeBiotSavart::REAL
            );
        A0Im.boundaryField()[control.interfaceLabel()] ==
            eBiotSavart.boundaryPatchA
            (
                control.interfaceLabel(),
                edgeBiotSavart::IMAG
            );
// TODO: Read nCorr from dict
        // Solve A0 in conductor region
        for (int corr=0; corr<5; corr++)
        {
            // Prepare A0-equations
            fvVectorMatrix A0ReEqn(fvm::laplacian(A0Re));
            fvVectorMatrix A0ImEqn(fvm::laplacian(A0Im));

            A0ReEqn.solve();
            A0ImEqn.solve();

            A0Re.correctBoundaryConditions();
            A0Im.correctBoundaryConditions();
        }

        manager.messages().newLine();

        // Reverse map A0 to default region
        globalStorage.A0Re().rmap(CONDUCTOR);
        globalStorage.A0Im().rmap(CONDUCTOR);

        eBiotSavart.writeVTK();
    }
