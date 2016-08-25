
    // Update sigma and A in conductor region
    // from data in default region
    {
        SM_REGIONSCOPE(CONDUCTOR);

        // Update sigma in conductor region
        globalStorage.sigma().mapCopyInternal(CONDUCTOR);

        // Update A in conductor region
        globalStorage.ARe().mapExtrapolate(CONDUCTOR);
        globalStorage.AIm().mapExtrapolate(CONDUCTOR);

        // Update rate of change of flux of A
        {
            storage.ddtARePhi() =  omega0*(fvc::interpolate(globalStorage.ARe()[CONDUCTOR]) & mesh.Sf());
            storage.ddtAImPhi() = -omega0*(fvc::interpolate(globalStorage.AIm()[CONDUCTOR]) & mesh.Sf());

            if (globalSettings.biotSavart)
            {
                storage.ddtARePhi() += omega0*(fvc::interpolate(globalStorage.A0Re()[CONDUCTOR]) & mesh.Sf());
                storage.ddtAImPhi() -= omega0*(fvc::interpolate(globalStorage.A0Im()[CONDUCTOR]) & mesh.Sf());
            }

            storage.ddtARePhi().correctBoundaryConditions();
            storage.ddtAImPhi().correctBoundaryConditions();
        }
    }
