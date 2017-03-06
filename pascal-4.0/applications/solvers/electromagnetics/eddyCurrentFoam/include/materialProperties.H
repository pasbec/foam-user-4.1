
    // Update sigma in default region
    {
        if (Control::debug)
        {
            Info<< Control::typeName << " | sigmaUpdate.H : "
                << "Update sigma."
                << endl;
        }

        // Correct boundary conditions of sigma
        // NOTE: Important for parallel runs
        globalStorage.sigma()[CONDUCTOR].correctBoundaryConditions();

        // Update sigma
        globalStorage.sigma().rmap(CONDUCTOR);
    }

    // Update mur in default region
    {
        if (Control::debug)
        {
            Info<< Control::typeName << " | murUpdate.H : "
                << "Update mur."
                << endl;
        }

        // Correct boundary conditions of mur
        // NOTE: Important for parallel runs
        globalStorage.mur()[CONDUCTOR].correctBoundaryConditions();

        // Update mur
        globalStorage.mur().rmap(CONDUCTOR);

        // Update nur
        globalStorage.nur()[CONDUCTOR] = 1.0/globalStorage.mur()[CONDUCTOR];
        globalStorage.nur().rmap(CONDUCTOR);
    }


