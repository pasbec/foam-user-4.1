
    // TODO: Remove interfaceLabel? Use alias of calculatedFvPatchField or special patchField

    // Solve V
    //
    // Hint: The electric scalar potential is only necessary
    // in case of a three-dimensional mesh
    //
    if (control.meshIs3D())
    {
        SM_REGIONSCOPE(CONDUCTOR);

        volScalarField& sigma = globalStorage.sigma()[CONDUCTOR];

        label complexSign = 0;
        volScalarField* VPtr = NULL;
        volVectorField* APtr = NULL;
        volVectorField* A0Ptr = NULL;
        surfaceScalarField* ddtAPhiPtr = NULL;

        if (control.switchV())
        {
            complexSign = 1;
            VPtr = &storage.VRe();
            APtr = &globalStorage.AIm()[CONDUCTOR];
            ddtAPhiPtr = &storage.ddtAImPhi();

            if (globalSettings.biotSavart)
            {
                A0Ptr = &globalStorage.A0Im()[CONDUCTOR];
            }
        }
        else
        {
            complexSign = -1;
            VPtr = &storage.VIm();
            APtr = &globalStorage.ARe()[CONDUCTOR];
            ddtAPhiPtr = &storage.ddtARePhi();

            if (globalSettings.biotSavart)
            {
                A0Ptr = &globalStorage.A0Re()[CONDUCTOR];
            }
        }

        volScalarField& V = *VPtr;
        volVectorField& A = *APtr;
        surfaceScalarField& ddtAPhi = *ddtAPhiPtr;

// TODO: Use dedicated boundary condition?
        // Update gradient boundary condition for V
        {
            label patchI = control.interfaceLabel();

            if
            (
                V.boundaryField()[patchI].type()
             == fixedGradientFvPatchField<scalar>::typeName
            )
            {
                fixedGradientFvPatchField<scalar>& VP =
                    refCast<fixedGradientFvPatchField<scalar> >
                    (
                        V.boundaryField()[patchI]
                    );

                // Set correct gradient
                VP.gradient() = -ddtAPhi.boundaryField()[patchI]/VP.patch().magSf();
            }

            V.correctBoundaryConditions();
        }

        // Prepare V-equation
        fvScalarMatrix VEqn(fvm::laplacian(sigma, V));

        // TODO: Apply electricPotentialLaplacian only for given faceSet?
        //       In this case we'd need a special electricConducivityGrad
        //       which evaluates to zero in normal direction to corresponding
        //       faces
        if
        (
            word
            (
                mesh.schemesDict().laplacianScheme
                (
                    "laplacian(" + sigma.name() + "," + V.name() + ")"
                )
            ) != "electricPotentialLaplacian"
        )
        {
            volVectorField sigmaGrad = fvc::grad(sigma);

            // Explicit source term for V due to variable sigma
            // in conducting region with implied condition div(A)=0
            VEqn -= complexSign * omega0 * (sigmaGrad & A);
//          VEqn -= complexSign * omega0 * sigma * fvc::div(A);

            if (globalSettings.biotSavart)
            {
                volVectorField& A0 = *A0Ptr;

                VEqn -= complexSign * omega0 * (sigmaGrad & A0);
            }
        }

        // Solve V-equation
        // NOTE: Currently a deflated solver is necessary
        //       as both matrices are kept singluar!
// TODO: Warning if no deflated solver is used!
        VEqn.solve(control.subDict(Region::CONDUCTOR, V.name()));

        // Correct boundary conditions
        V.correctBoundaryConditions();
    }
