
    // Update flux of gradient of normal component of A
    // in default region
    {
        SM_REGIONSCOPE(DEFAULT);

        for (int cI = 0; cI < 2; cI++)
        {
            volVectorField* APtr = NULL;
            surfaceVectorField* phiGradAnPtr = NULL;

            if (cI == 0)
            {
                APtr = &globalStorage.ARe()[DEFAULT];
                phiGradAnPtr = &storage.phiGradAnRe();
            }
            else
            {
                APtr = &globalStorage.AIm()[DEFAULT];
                phiGradAnPtr = &storage.phiGradAnIm();
            }

            volVectorField& A = *APtr;
            surfaceVectorField& phiGradAn = *phiGradAnPtr;

            surfaceTensorField gradAf =
                fvc::interpolate
                (
                    fvc::grad(linearInterpolate(A)),
                    "interpolate(grad(" + A.name() + "))"
                );

            phiGradAn = -gradAf & mesh.Sf();

            phiGradAn += tr(gradAf) * mesh.Sf();
        }
    }

    // Update A in conductor region
    // from data in default region
    {
        SM_REGIONSCOPE(CONDUCTOR);

        if (Control::debug)
        {
            Info<< Control::typeName << " | AUpdate.H : "
                << "Update A in conductor region."
                << endl;
        }

// TODO: Use relaxation?
//         globalStorage.ARe()[CONDUCTOR].storePrevIter();
//         globalStorage.ARe()[CONDUCTOR].storePrevIter();

        // Update A
        globalStorage.ARe().mapInterpolate(CONDUCTOR);
        globalStorage.AIm().mapInterpolate(CONDUCTOR);

// TODO: Use relaxation?
//         const scalar AReRf = globalMesh[DEFAULT].solutionDict().fieldRelaxationFactor(globalStorage.ARe()[DEFAULT].name());
//         const scalar AImRf = globalMesh[DEFAULT].solutionDict().fieldRelaxationFactor(globalStorage.AIm()[DEFAULT].name());
//         globalStorage.ARe()[CONDUCTOR].relax(AReRf);
//         globalStorage.ARe()[CONDUCTOR].relax(AImRf);
    }


    // Update flux of rate of change of A
    // in conductor region
    if (control.meshIs3D())
    {
        SM_REGIONSCOPE(CONDUCTOR);

// TEST
// TODO: Interpolation with nuf?
        storage.phiDdtARe() =  omega0*(fvc::interpolate(globalStorage.ARe()[CONDUCTOR]) & mesh.Sf());
        storage.phiDdtAIm() = -omega0*(fvc::interpolate(globalStorage.AIm()[CONDUCTOR]) & mesh.Sf());

        if (globalSettings.biotSavart)
        {
            storage.phiDdtARe() += omega0*(fvc::interpolate(globalStorage.A0Re()[CONDUCTOR]) & mesh.Sf());
            storage.phiDdtAIm() -= omega0*(fvc::interpolate(globalStorage.A0Im()[CONDUCTOR]) & mesh.Sf());
        }
    }
