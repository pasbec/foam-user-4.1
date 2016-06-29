
    // Solve V
    //
    // Hint: The electric scalar potential is only necessary
    // in case of a three-dimensional mesh
    //
    if (control.meshIs3D())
    {
        SM_REGIONSCOPE(CONDUCTOR);

        volScalarField& sigma = globalStorage.sigma()[CONDUCTOR];

        volVectorField& ARe = globalStorage.ARe()[CONDUCTOR];
        volVectorField& AIm = globalStorage.AIm()[CONDUCTOR];

        surfaceScalarField& ddtARePhi = storage.ddtARePhi();
        surfaceScalarField& ddtAImPhi = storage.ddtAImPhi();

        volScalarField& VRe = storage.VRe();
        volScalarField& VIm = storage.VIm();

// TODO: Use dedicated boundary condition?
        // Update gradient boundary condition for V
        {
            label patchI = control.interfaceLabel();

            if
            (
                VRe.boundaryField()[patchI].type()
             == fixedGradientFvPatchField<scalar>::typeName
            )
            {
                fixedGradientFvPatchField<scalar>& VReP =
                    refCast<fixedGradientFvPatchField<scalar> >
                    (
                        VRe.boundaryField()[patchI]
                    );

                // Set correct gradient
                VReP.gradient() = -ddtAImPhi.boundaryField()[patchI]/VReP.patch().magSf();
            }

            if
            (
                VIm.boundaryField()[patchI].type()
             == fixedGradientFvPatchField<scalar>::typeName
            )
            {
                fixedGradientFvPatchField<scalar>& VImP =
                    refCast<fixedGradientFvPatchField<scalar> >
                    (
                        VIm.boundaryField()[patchI]
                    );

                // Set correct gradient
                VImP.gradient() = -ddtARePhi.boundaryField()[patchI]/VImP.patch().magSf();
            }

            VRe.correctBoundaryConditions();
            VIm.correctBoundaryConditions();
        }

        // Prepare V-equations
        fvScalarMatrix VReEqn(fvm::laplacian(sigma, VRe));
        fvScalarMatrix VImEqn(fvm::laplacian(sigma, VIm));

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
                    "laplacian(" + sigma.name() + "," + VRe.name() + ")"
                )
            ) != "electricPotentialLaplacian"
        )
        {
            volVectorField sigmaGrad = fvc::grad(sigma);

            // Explicit source term for VRe due to variable sigma
            // in conducting region with implied condition div(A)=0
            VReEqn -= omega0 * (sigmaGrad & AIm);
//          VReEqn -= omega0 * sigma * fvc::div(AIm);

            if (globalSettings.biotSavart)
            {
                VReEqn -= omega0 * (sigmaGrad & globalStorage.A0Im()[CONDUCTOR]);
//              VReEqn -= omega0 * sigma * fvc::div(globalStorage.A0Im()[CONDUCTOR]);
            }
        }

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
                    "laplacian(" + sigma.name() + "," + VIm.name() + ")"
                )
            ) != "electricPotentialLaplacian"
        )
        {
            volVectorField sigmaGrad = fvc::grad(sigma);

            // Explicit source term for VIm due to variable sigma
            // in conducting region with implied condition div(A)=0
            VImEqn += omega0 * (sigmaGrad & ARe);
//          VImEqn += omega0 * sigma * fvc::div(ARe);

            if (globalSettings.biotSavart)
            {
                VImEqn += omega0 * (sigmaGrad & globalStorage.A0Re()[CONDUCTOR]);
//              VImEqn += omega0 * sigma * fvc::div(globalStorage.A0Re()[CONDUCTOR]);
            }
        }

        // Solve V-equations
        // NOTE: Currently a deflated solver is necessary
        //       as both matrices are kept singluar!
// TODO: Warning if no deflated solver is used!
        VReEqn.solve(control.subDict(Region::CONDUCTOR, VRe.name()));
        VImEqn.solve(control.subDict(Region::CONDUCTOR, VIm.name()));

        // Correct boundary conditions
        VRe.correctBoundaryConditions();
        VIm.correctBoundaryConditions();
    }
