
    // Init gradient of V in first iteration
    if (control.meshIs3D() )
    {
        SM_REGIONSCOPE(CONDUCTOR);

        globalStorage.VReGrad()[CONDUCTOR] == fvc::grad(storage.VRe());
        globalStorage.VImGrad()[CONDUCTOR] == fvc::grad(storage.VIm());

        globalStorage.VReGrad().rmap(CONDUCTOR);
        globalStorage.VImGrad().rmap(CONDUCTOR);
    }

    // AV iterations
    while (control.loop())
    {
        // Solve for A in base region
        {
            SM_REGIONSCOPE(DEFAULT);

            volScalarField& sigma = globalStorage.sigma()[DEFAULT];

            volVectorField& ARe = globalStorage.ARe()[DEFAULT];
            volVectorField& AIm = globalStorage.AIm()[DEFAULT];

            volVectorField& VReGrad = globalStorage.VReGrad()[DEFAULT];
            volVectorField& VImGrad = globalStorage.VImGrad()[DEFAULT];

            volVectorField& j0Re = storage.j0Re();
            volVectorField& j0Im = storage.j0Im();

#           include "AEqn.H"
        }

        // The electric scalar potential is only necesasry if have
        // a 3D case
        if (control.meshIs3D())
        {
            // Update sigma in conductor region
            globalStorage.sigma().mapExtrapolate(CONDUCTOR);

            // Update A in conductor region
            globalStorage.ARe().mapExtrapolate(CONDUCTOR);
            globalStorage.AIm().mapExtrapolate(CONDUCTOR);

            // Solve for V in conductor region
            {
                SM_REGIONSCOPE(CONDUCTOR);
// TODO
                // Relax A
//                 control.relax(globalStorage.ARe());
//                 control.relax(globalStorage.AIm());

                volScalarField& sigma = globalStorage.sigma()[CONDUCTOR];

                volVectorField& ARe = globalStorage.ARe()[CONDUCTOR];
                volVectorField& AIm = globalStorage.AIm()[CONDUCTOR];

                volScalarField& VRe = storage.VRe();
                volScalarField& VIm = storage.VIm();

#               include "VEqn.H"
            }

            // Update gradient of V
            {
                SM_REGIONSCOPE(CONDUCTOR);

                globalStorage.VReGrad()[CONDUCTOR] == fvc::grad(storage.VRe());
                globalStorage.VImGrad()[CONDUCTOR] == fvc::grad(storage.VIm());

                globalStorage.VReGrad().rmap(CONDUCTOR);
                globalStorage.VImGrad().rmap(CONDUCTOR);
            }
        }

        manager.messages().newLine();
    }
