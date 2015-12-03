
    // Solve V
    {
        // Update gradient boundary condition for V
        {
            const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();

            label patchI = boundaryMesh.findPatchID(control.interfaceName());

            if (patchI == -1)
            {
                FatalError
                    << "Given patch name " << control.interfaceName()
                    << " for point mapping does not exist."
                    << abort(FatalError);
            }

            polyPatch patch = mesh.boundaryMesh()[patchI];

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
                    VReP.gradient() =  control.omega().value()
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
                    VImP.gradient() = -control.omega().value()
                    * (ARe.boundaryField()[patchI] & VImP.patch().nf());

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
        VReEqn -= control.omega() * (sigmaGrad & AIm);
        VImEqn += control.omega() * (sigmaGrad & ARe);

//         // TEST
//         VReEqn -= control.omega() * sigma * fvc::div(AIm);
//         VImEqn += control.omega() * sigma * fvc::div(ARe);

//         // Restore solution level
//         cVRe += VReAverage;
//         cVIm += VImAverage;

        // Solve V-equations
        // NOTE: Currently a deflated solver is necessary
        // TODO: WARNING if no deflated solver is used!
        //       as both matrices are kept singluar!
        lduSolverPerformance VReEqnPerf = VReEqn.solve(control.Vdict());
        lduSolverPerformance VImEqnPerf = VImEqn.solve(control.Vdict());

        // Store max residual
        scalar resVRe = VReEqnPerf.initialResidual();
        scalar resVIm = VImEqnPerf.initialResidual();
        control.setResidualOfV(max(resVRe, resVIm));

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
