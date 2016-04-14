
    // Solve V
    {
        SM_REGIONSCOPE(CONDUCTOR);

        volScalarField& sigma = globalStorage.sigma()[CONDUCTOR];

        volVectorField& ARe = globalStorage.ARe()[CONDUCTOR];
        volVectorField& AIm = globalStorage.AIm()[CONDUCTOR];

        volScalarField& VRe = storage.VRe();
        volScalarField& VIm = storage.VIm();

        // Update gradient boundary condition for V
        {
            const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();

            label patchI = boundaryMesh.findPatchID(control.interfaceName());

            if (patchI == -1)
            {
                FatalError
                    << "Given patch name "
                    << control.interfaceName()
                    << " does not exist."
                    << abort(FatalError);
            }

            polyPatch patch = mesh.boundaryMesh()[patchI];

// TODO: Use seperate boundary condition
            // Boundary values for V
            {
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
                    VReP.gradient() =  omega0.value()
                      * (AIm.boundaryField()[patchI] & VReP.patch().nf());

                    // Add part of A0 if Biot-Savert is used
                    if (globalSettings.biotSavart)
                    {
                        VReP.gradient() += omega0.value()
                          * (globalStorage.A0Im()[CONDUCTOR].boundaryField()[patchI]
                          & VReP.patch().nf());
                    }

                    // Evaluate values
                    VReP.evaluate();
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
                    VImP.gradient() = -omega0.value()
                      * (ARe.boundaryField()[patchI] & VImP.patch().nf());

                    // Add part of A0 if Biot-Savert is used
                    if (globalSettings.biotSavart)
                    {
                        VImP.gradient() -= omega0.value()
                          * (globalStorage.A0Re()[CONDUCTOR].boundaryField()[patchI]
                          & VImP.patch().nf());
                    }

                    // Evaluate values
                    VImP.evaluate();
                }
            }

            // Correct boundary conditions for V
            VRe.correctBoundaryConditions();
            VIm.correctBoundaryConditions();
        }

        // Prepare V-equations
        fvScalarMatrix VReEqn(fvm::laplacian(sigma, VRe));
        fvScalarMatrix VImEqn(fvm::laplacian(sigma, VIm));

        volVectorField sigmaGrad = fvc::grad(sigma);

        // Explicit source terms due to variable sigma
        // in conducting region with implied condition div(A)=0
        VReEqn -= omega0 * (sigmaGrad & AIm);
        VImEqn += omega0 * (sigmaGrad & ARe);

//         // TEST
//         VReEqn -= omega0 * sigma * fvc::div(AIm);
//         VImEqn += omega0 * sigma * fvc::div(ARe);

        if (globalSettings.biotSavart)
        {
            VReEqn -= omega0 * (sigmaGrad & globalStorage.A0Im()[CONDUCTOR]);
            VImEqn += omega0 * (sigmaGrad & globalStorage.A0Re()[CONDUCTOR]);

//             // TEST
//             VReEqn -= omega0 * sigma * fvc::div(globalStorage.A0Im()[CONDUCTOR]);
//             VImEqn += omega0 * sigma * fvc::div(globalStorage.A0Re()[CONDUCTOR]);
        }

//         // Restore solution level from average
//         cVRe += VReAverage;
//         cVIm += VImAverage;

        // Solve V-equations
        // NOTE: Currently a deflated solver is necessary
        //       as both matrices are kept singluar!
// TODO: Warning if no deflated solver is used!
        VReEqn.solve(control.subDict(Region::CONDUCTOR, VRe.name()));
        VImEqn.solve(control.subDict(Region::CONDUCTOR, VIm.name()));

        // Correct boundary conditions
        VRe.correctBoundaryConditions();
        VIm.correctBoundaryConditions();

//         // Calculate average
//         VReAverage = average(cVRe);
//         VImAverage = average(cVIm);
//
//         // Normalize solution level
//         cVRe -= VReAverage;
//         cVIm -= VImAverage;
    }
