
    Info << "Total surface tension force: "
         << interface.totalSurfaceTensionForce() << endl;

    vector totalForce =
        interface.totalViscousForce()
        + interface.totalPressureForce();

    Info << "Total force: " << totalForce << nl << endl;