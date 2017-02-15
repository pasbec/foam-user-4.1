
    // Update flux of gradient of normal component of A0
    // in default region
    if (globalSettings.biotSavart)
    {
        SM_REGIONSCOPE(CONDUCTOR);

        for (int cI = 0; cI < 2; cI++)
        {
            volVectorField* A0Ptr = NULL;
            surfaceVectorField* phiGradA0nPtr = NULL;

            if (cI == 0)
            {
                A0Ptr = &globalStorage.A0Re()[CONDUCTOR];
                phiGradA0nPtr = &storage.phiGradA0nRe();
            }
            else
            {
                A0Ptr = &globalStorage.A0Im()[CONDUCTOR];
                phiGradA0nPtr = &storage.phiGradA0nIm();
            }

            volVectorField& A0 = *A0Ptr;
            surfaceVectorField& phiGradA0n = *phiGradA0nPtr;

            surfaceTensorField gradAf =
                fvc::interpolate
                (
                    fvc::grad(linearInterpolate(A0)),
                    "interpolate(grad(" + A0.name() + "))"
                );

            phiGradA0n = -gradAf & mesh.Sf();

            phiGradA0n += tr(gradAf) * mesh.Sf();
        }
    }

    // Update A0 in default region
    // from data in conducting region
    if (globalSettings.biotSavart)
    {
        SM_REGIONSCOPE(DEFAULT);

        if (Control::debug)
        {
            Info<< Control::typeName << " | AUpdate.H : "
                << "Update A0 in default region."
                << endl;
        }

        // Update A0
        globalStorage.A0Re().rmap(CONDUCTOR);
        globalStorage.A0Im().rmap(CONDUCTOR);
    }
