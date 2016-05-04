
    // Update sigma and A in conductor region
    // from data in default region
    {
        SM_REGIONSCOPE(CONDUCTOR);

        globalStorage.sigma().mapExtrapolate(CONDUCTOR);

        // Update A in conductor region
// TODO: Extrapolation
//         globalStorage.ARe().mapExtrapolate(CONDUCTOR);
//         globalStorage.AIm().mapExtrapolate(CONDUCTOR);
        globalStorage.ARe().mapCopyInternal(CONDUCTOR);
        globalStorage.AIm().mapCopyInternal(CONDUCTOR);

        // Update rate of change of flux of A
        {
            storage.ddtARePhi() =  omega0*(fvc::interpolate(globalStorage.ARe()[CONDUCTOR]) & mesh.Sf());
            storage.ddtAImPhi() = -omega0*(fvc::interpolate(globalStorage.AIm()[CONDUCTOR]) & mesh.Sf());

            if (globalSettings.biotSavart)
            {
                storage.ddtARePhi() += omega0*(fvc::interpolate(globalStorage.A0Re()[CONDUCTOR]) & mesh.Sf());
                storage.ddtAImPhi() -= omega0*(fvc::interpolate(globalStorage.A0Im()[CONDUCTOR]) & mesh.Sf());
            }
        }
    }
