
    // Update gradient of V in default region
    // from data in conductor region
    if (control.meshIs3D())
    {
        if (Control::debug)
        {
            Info<< Control::typeName << " | BUpdate.H : "
                << "Update gradient of V."
                << endl;
        }

        SM_REGIONSCOPE(CONDUCTOR);

        volVectorField& VReGrad = globalStorage.VReGrad()[CONDUCTOR];
        volVectorField& VImGrad = globalStorage.VImGrad()[CONDUCTOR];

// TODO: Use relaxation?
//         VReGrad.storePrevIter();
//         VImGrad.storePrevIter();

        // Calculate gradient of V
        VReGrad = fvc::grad(storage.VRe());
        VImGrad = fvc::grad(storage.VIm());

// TODO: Use relaxation?
//         VReGrad.relax(0.5);
//         VImGrad.relax(0.5);

        // Update gradient of V
        globalStorage.VReGrad().rmap(CONDUCTOR);
        globalStorage.VImGrad().rmap(CONDUCTOR);
    }
