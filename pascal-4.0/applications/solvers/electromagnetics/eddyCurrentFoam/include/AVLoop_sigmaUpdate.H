
    // Update sigma in default region
    {
        SM_REGIONSCOPE(CONDUCTOR);

        if (Control::debug)
        {
            Info<< Control::typeName << " | sigmaUpdate.H : "
                << "Update sigma."
                << endl;
        }

        // Update sigma
        globalStorage.sigma().rmap(CONDUCTOR);
    }
