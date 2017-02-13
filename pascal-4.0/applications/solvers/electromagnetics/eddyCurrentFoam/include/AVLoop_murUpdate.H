
    // Update mur in default region
    {
        SM_REGIONSCOPE(CONDUCTOR);

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

        // Update rMur
        globalStorage.rMur()[CONDUCTOR] = 1.0/globalStorage.mur()[CONDUCTOR];
        globalStorage.rMur().rmap(CONDUCTOR);
    }
