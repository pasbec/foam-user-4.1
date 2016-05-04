
    // Solve V
    {
        SM_REGIONSCOPE(CONDUCTOR);

        volScalarField& sigma = globalStorage.sigma()[CONDUCTOR];

        volVectorField& ARe = globalStorage.ARe()[CONDUCTOR];
        volVectorField& AIm = globalStorage.AIm()[CONDUCTOR];

        surfaceScalarField& ddtARePhi = storage.ddtARePhi();
        surfaceScalarField& ddtAImPhi = storage.ddtAImPhi();

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

// TODO: Use dedicated boundary condition
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
            }

            // Correct boundary conditions for V
            VRe.correctBoundaryConditions();
            VIm.correctBoundaryConditions();
        }

        // Prepare V-equations
        fvScalarMatrix VReEqn(fvm::laplacian(sigma, VRe));
        fvScalarMatrix VImEqn(fvm::laplacian(sigma, VIm));

// TODO FIXME: Think wether this is necessary or not...
//         volVectorField sigmaGrad = fvc::grad(sigma);
//
//         // Explicit source terms due to variable sigma
//         // in conducting region with implied condition div(A)=0
//         VReEqn -= omega0 * (sigmaGrad & AIm);
//         VImEqn += omega0 * (sigmaGrad & ARe);
//
// //         // TEST
// //         VReEqn -= omega0 * sigma * fvc::div(AIm);
// //         VImEqn += omega0 * sigma * fvc::div(ARe);
//
//         if (globalSettings.biotSavart)
//         {
//             VReEqn -= omega0 * (sigmaGrad & globalStorage.A0Im()[CONDUCTOR]);
//             VImEqn += omega0 * (sigmaGrad & globalStorage.A0Re()[CONDUCTOR]);
//
// //             // TEST
// //             VReEqn -= omega0 * sigma * fvc::div(globalStorage.A0Im()[CONDUCTOR]);
// //             VImEqn += omega0 * sigma * fvc::div(globalStorage.A0Re()[CONDUCTOR]);
//         }

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
