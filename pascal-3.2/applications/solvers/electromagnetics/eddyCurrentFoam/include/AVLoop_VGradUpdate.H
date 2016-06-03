
    // Update gradient of V in default region
    // from data in conductor region
    {
        SM_REGIONSCOPE(CONDUCTOR);

        globalStorage.VReGrad()[CONDUCTOR] = fvc::grad(storage.VRe());
        globalStorage.VImGrad()[CONDUCTOR] = fvc::grad(storage.VIm());

        globalStorage.VReGrad().rmap(CONDUCTOR);
        globalStorage.VImGrad().rmap(CONDUCTOR);
    }
