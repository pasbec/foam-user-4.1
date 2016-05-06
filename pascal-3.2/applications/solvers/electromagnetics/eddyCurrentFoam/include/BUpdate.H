
    // Calculate B
    {
        globalStorage.BRe()[DEFAULT] == fvc::curl(globalStorage.ARe()[DEFAULT]);
        globalStorage.BIm()[DEFAULT] == fvc::curl(globalStorage.AIm()[DEFAULT]);

// TODO: Extrapolation
//         globalStorage.BRe().mapExtrapolate(CONDUCTOR);
//         globalStorage.BIm().mapExtrapolate(CONDUCTOR);
        globalStorage.BRe().mapCopyInternal(CONDUCTOR);
        globalStorage.BIm().mapCopyInternal(CONDUCTOR);

        // If Biot-Savart IS used, A0 in the non-cunducting region
        // is NOT calculated. Thus, we are not allowed to do fvc::curl(A + A0)
        // in the whole default region.
        if (globalSettings.biotSavart)
        {
            globalStorage.BRe()[CONDUCTOR] += fvc::curl(globalStorage.A0Re()[CONDUCTOR]);
            globalStorage.BIm()[CONDUCTOR] += fvc::curl(globalStorage.A0Im()[CONDUCTOR]);

            globalStorage.BRe().rmap(CONDUCTOR);
            globalStorage.BIm().rmap(CONDUCTOR);
        }
    }
