
    // Update sigma and A in conductor region
    // from data in default region
    {
        SM_REGIONSCOPE(CONDUCTOR);

// TODO: Extrapolation bugfix for skewed meshes
        globalStorage.sigma().mapExtrapolate(CONDUCTOR);
//         globalStorage.sigma().mapCopyInternal(CONDUCTOR);

        // Update A in conductor region
        if (globalSettings.biotSavart)
        {
// TODO: Extrapolation bugfix for skewed meshes
            globalStorage.ARe().mapExtrapolate(CONDUCTOR);
            globalStorage.AIm().mapExtrapolate(CONDUCTOR);
//             globalStorage.ARe().mapCopyInternal(CONDUCTOR);
//             globalStorage.AIm().mapCopyInternal(CONDUCTOR);
        }
        else
        {
// TODO: Interpolation in parallel
            globalStorage.ARe().mapInterpolate(CONDUCTOR);
            globalStorage.AIm().mapInterpolate(CONDUCTOR);
//             globalStorage.ARe().mapCopyInternal(CONDUCTOR);
//             globalStorage.AIm().mapCopyInternal(CONDUCTOR);
        }

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
