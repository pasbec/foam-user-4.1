
    // Update sigma in default region
    {
        SM_REGIONSCOPE(CONDUCTOR);

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
