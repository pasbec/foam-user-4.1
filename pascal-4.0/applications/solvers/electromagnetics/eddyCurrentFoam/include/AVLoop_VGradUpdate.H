
    // Update gradient of V in default region
    // from data in conductor region
    //
    // Hint: The electric scalar potential is only necessary
    // in case of a three-dimensional mesh
    //
    if (control.meshIs3D())
    {
        SM_REGIONSCOPE(CONDUCTOR);

        globalStorage.VReGrad()[CONDUCTOR] = fvc::grad(storage.VRe());
        globalStorage.VImGrad()[CONDUCTOR] = fvc::grad(storage.VIm());

        globalStorage.VReGrad().rmap(CONDUCTOR);
        globalStorage.VImGrad().rmap(CONDUCTOR);
    }
