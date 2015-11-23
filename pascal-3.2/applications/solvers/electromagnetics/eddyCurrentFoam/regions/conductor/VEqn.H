
    // Solve V
    {
        // Update gradient boundary condition for V
        {
            forAll(mesh.boundaryMesh(), patchI)
            {
                polyPatch patch = mesh.boundaryMesh()[patchI];

                // Old internal faces
                // TODO: What about an more appropriate patch name?
                if (patch.name() == "oldInternalFaces")
                {
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
                            VReP.gradient() =  eddy.omega().value()
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
                            VImP.gradient() = -eddy.omega().value()
                            * (ARe.boundaryField()[patchI] & VImP.patch().nf());

                            // Evaluate values
                            VImP.evaluate();
                        }
                    }
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
        VReEqn -= eddy.omega() * (sigmaGrad & AIm);
        VImEqn += eddy.omega() * (sigmaGrad & ARe);

//         // TEST
//         VReEqn -= eddy.omega() * sigma * fvc::div(AIm);
//         VImEqn += eddy.omega() * sigma * fvc::div(ARe);

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
    }
