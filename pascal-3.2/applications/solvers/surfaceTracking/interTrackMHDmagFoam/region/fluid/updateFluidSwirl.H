
    // Remove swirl
    if (removeSwirl > -1)
    {
        Info << "Removing swirl-velocity-component: " << removeSwirl << nl;
        U.field().replace(removeSwirl, 0.0);
        U.boundaryField().replace(removeSwirl, 0.0);
    }