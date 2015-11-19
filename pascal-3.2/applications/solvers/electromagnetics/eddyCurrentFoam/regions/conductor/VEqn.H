
    // Solve V
    {

        // Boundary fields
        forAll(mesh.boundaryMesh(), patchI)
        {
            polyPatch patch = mesh.boundaryMesh()[patchI];

// FIXME: START
            // Create interpolated surface field of A
            surfaceVectorField sARe("sARe", fvc::interpolate(ARe_[defaultRegionID]));
            surfaceVectorField sAIm("sAIm", fvc::interpolate(AIm_[defaultRegionID]));

            // Get ref to face map
            labelList meshFaceMap = mesh_.faceMap(conductorRegionID);
// FIXME: END

            // Old internal faces
            // TODO: What about an more appropriate patch name?
            if (patch.name() == "oldInternalFaces")
            {
// FIXME: START

                // TODO: ARe/AIm
                // Boundary values for A
                {
                    vectorField& AReP = ARe.boundaryField()[patchI];
                    vectorField& AImP = AIm.boundaryField()[patchI];

                    forAll(patch, facei)
                    {
                        label faceI = patch.start() + facei;

                        AReP[facei] = sARe[meshFaceMap[faceI]];
                        AImP[facei] = sAIm[meshFaceMap[faceI]];
                    }
                }

                // TODO: VRe/VIm
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
                        VReP.gradient() =  omega.value()
                        * (AIm.boundaryField()[patchI] & VReP.patch().nf());

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
                        VImP.gradient() = -omega.value()
                        * (ARe.boundaryField()[patchI] & VImP.patch().nf());

                        // Evaluate values
                        VImP.evaluate();
                    }
                }
// FIXME: END

                // Boundary values for sigma
                {
                    fvPatchField<scalar>& sigmaP = sigma.boundaryField()[patchI];

                    // Evaluate values
                    sigmaP = sigmaP.patchInternalField();
                }
            }

            // TODO: What about other patches?
        }

        // Correct boundary conditions for A
        ARe.correctBoundaryConditions();
        AIm.correctBoundaryConditions();

        // Correct boundary conditions for V
        VRe.correctBoundaryConditions();
        VIm.correctBoundaryConditions();

        // Correct boundary conditions for sigma
        sigma.correctBoundaryConditions();

        // Prepare V-equations
        fvScalarMatrix VReEqn(fvm::laplacian(sigma, VRe));
        fvScalarMatrix VImEqn(fvm::laplacian(sigma, VIm));

        volVectorField sigmaGrad = fvc::grad(sigma);

        // Explicit source terms due to variable sigma
        // in conducting region with implied condition div(A)=0
        VReEqn -= omega * (sigmaGrad & AIm);
        VImEqn += omega * (sigmaGrad & ARe);

//         // TEST
//         VReEqn -= omega * sigma * fvc::div(AIm);
//         VImEqn += omega * sigma * fvc::div(ARe);

//         // Restore solution level
//         cVRe += VReAverage;
//         cVIm += VImAverage;

        // Solve V-equations
        // NOTE: Currently a deflated solver is necessary
        //       as both matrices are kept singluar!
        lduSolverPerformance VReEqnPerf = VReEqn.solve(Vdict);
        lduSolverPerformance VImEqnPerf = VImEqn.solve(Vdict);

        // Store max residual
        scalar resVRe = VReEqnPerf.initialResidual();
        scalar resVIm = VImEqnPerf.initialResidual();
        resV = max(resVRe, resVIm);

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

        // Calculate gradient
        VReGrad = fvc::grad(VRe);
        VImGrad = fvc::grad(VIm);
    }
