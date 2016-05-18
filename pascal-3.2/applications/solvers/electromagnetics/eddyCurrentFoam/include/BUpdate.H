
    // Calculate B
    {
        globalStorage.BRe()[DEFAULT] == fvc::curl(globalStorage.ARe()[DEFAULT]);
        globalStorage.BIm()[DEFAULT] == fvc::curl(globalStorage.AIm()[DEFAULT]);

        // Update B in conductor region
        if (globalSettings.biotSavart)
        {
// TODO: Extrapolation bugfix for skewed meshes
            globalStorage.BRe().mapExtrapolate(CONDUCTOR);
            globalStorage.BIm().mapExtrapolate(CONDUCTOR);
//             globalStorage.BRe().mapCopyInternal(CONDUCTOR);
//             globalStorage.BIm().mapCopyInternal(CONDUCTOR);
        }
        else
        {
// TODO: Interpolation in parallel
            globalStorage.BRe().mapInterpolate(CONDUCTOR);
            globalStorage.BIm().mapInterpolate(CONDUCTOR);
//             globalStorage.BRe().mapCopyInternal(CONDUCTOR);
//             globalStorage.BIm().mapCopyInternal(CONDUCTOR);
        }

        // If Biot-Savart IS used, A0 in the non-cunducting region
        // is NOT calculated. Thus, we are not allowed to do fvc::curl(A + A0)
        // in the whole default region.
        if (globalSettings.biotSavart)
        {
            globalStorage.BRe()[CONDUCTOR] += fvc::curl(globalStorage.A0Re()[CONDUCTOR]);
            globalStorage.BIm()[CONDUCTOR] += fvc::curl(globalStorage.A0Im()[CONDUCTOR]);

            globalStorage.BRe()[DEFAULT] *= 0;
            globalStorage.BIm()[DEFAULT] *= 0;

            globalStorage.BRe().rmap(CONDUCTOR);
            globalStorage.BIm().rmap(CONDUCTOR);
        }
    }
