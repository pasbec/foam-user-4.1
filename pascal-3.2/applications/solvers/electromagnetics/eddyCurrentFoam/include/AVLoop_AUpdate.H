
    // Update sigma and A in conductor region
    // from data in default region
    {
        globalStorage.sigma().mapExtrapolate(CONDUCTOR);

        // Update A in conductor region
// TODO: Extrapolation
//         globalStorage.ARe().mapExtrapolate(CONDUCTOR);
//         globalStorage.AIm().mapExtrapolate(CONDUCTOR);
        globalStorage.ARe().mapCopyInternal(CONDUCTOR);
        globalStorage.AIm().mapCopyInternal(CONDUCTOR);
    }
