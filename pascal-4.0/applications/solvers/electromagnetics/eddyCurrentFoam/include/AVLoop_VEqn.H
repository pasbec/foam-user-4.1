
    // TODO: Remove interfaceLabel? Use alias of calculatedFvPatchField or special patchField

    // Solve V
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

            if (isA<fixedGradientFvPatchScalarField>(V.boundaryField()[patchI]))
            {
                fixedGradientFvPatchScalarField& VP =
                    refCast<fixedGradientFvPatchScalarField>
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


// // TODO: Check for fixesValue
//         // Compatibility check
//         {
//             bool singular = true;
//
//             forAll(mesh.boundary(), patchI)
//             {
//                 if
//                 (
//                     !V.boundaryField()[patchI].coupled()
//                  && V.boundaryField()[patchI].fixesValue()
//                 )
//                 {
//                     singular = false;
//                 }
//             }
//
//             if (singular)
//             {
//                 scalar sourceVolSum = sum(VEqn.source() * mesh.V());
//                 scalar fluxAreaSum = 0;
//                 forAll(mesh.boundary(), patchI)
//                 {
//                     if (isA<fixedGradientFvPatchScalarField>(V.boundaryField()[patchI]))
//                     {
//                         fixedGradientFvPatchScalarField& VP =
//                             refCast<fixedGradientFvPatchScalarField>
//                             (
//                                 V.boundaryField()[patchI]
//                             );
// // TODO: Calc fluxAreaSum
//                     }
//                 }
//
//                 if (Control::debug)
//                 {
//                     Pout<< "eddyCurrentApp::Control | VEqn.H : "
//                         << "Local compatibility check for " << V.name()
//                         << "Eqn  (domain/boundary): " << sourceVolSum
//                         << " / " << fluxAreaSum
//                         << endl;
//                 }
//
//                 scalar sourceVolSumGlobal = returnReduce(sourceVolSum, sumOp<scalar>());
//                 scalar fluxAreaSumGlobal = returnReduce(fluxAreaSum, sumOp<scalar>());
//
//                 Info<< "Global compatibility check for " << V.name()
//                     << "Eqn (domain/boundary): " << sourceVolSumGlobal
//                     << " / " << fluxAreaSumGlobal
//                     << endl;
//             }
//         }


        // Solve V-equation
        // NOTE: Currently a deflated solver is necessary
        //       as both matrices are kept singluar!
// TODO: Warning if no deflated solver is used!
        VEqn.solve(control.subDict(Region::CONDUCTOR, V.name()));

        // Correct boundary conditions
        V.correctBoundaryConditions();
    }
