
    // Solve A component-wise
    {
        // Assemble base-matrices for A
        fvVectorMatrix AReEqn
        (
            rMu0 * fvm::laplacian(ARe)
//        + alpha * AIm
//        - sigma * VReGrad
          + jsRe
        );
        fvVectorMatrix AImEqn
        (
            rMu0 * fvm::laplacian(AIm)
//        - alpha * ARe
//        - sigma * VImGrad
          + jsIm
        );

        // Create list with boundary types
        wordList ARePatchTypes(mesh.boundary().size());
        wordList AImPatchTypes(mesh.boundary().size());

        forAll(ARe.boundaryField(), patchI)
        {
            ARePatchTypes[patchI] =
                ARe.boundaryField()[patchI].type();
        }
        forAll(AIm.boundaryField(), patchI)
        {
            AImPatchTypes[patchI] =
                AIm.boundaryField()[patchI].type();
        }

        // Relax gradient of V
        {
            // Store old gradient
            VReGrad.storePrevIter();
            VImGrad.storePrevIter();

            // Get relaxation factor
            scalar relax = AVdict.lookupOrDefault<scalar>("relax",1.0);

            if (relax < (1.0-SMALL))
            {
                // Relaxation of gradient of V
                VReGrad.relax(relax);
                VImGrad.relax(relax);
            }
        }

        // Component residual
        scalar resAcmpt = VSMALL;

        // Loop over components
        for (direction cmpt = 0; cmpt < vector::nComponents; cmpt++)
        {
            // Create component fields
            volScalarField AReCmpt
            (
                IOobject
                (
                    "A" + word(vector::componentNames[cmpt]) + "Re",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                ARe.component(cmpt),
                ARePatchTypes
            );
            volScalarField AImCmpt
            (
                IOobject
                (
                    "A" + word(vector::componentNames[cmpt]) + "Im",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                AIm.component(cmpt),
                AImPatchTypes
            );

            // Correct/Evaluate boundary values
            AReCmpt.correctBoundaryConditions();
            AImCmpt.correctBoundaryConditions();

            // Create complex vector2 field
            volVector2Field ACmpt
            (
                IOobject
                (
                    "A" + word(vector::componentNames[cmpt]),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedVector2(word(),dimless,vector2::zero)
            );

            // Init complex vector2 field
            ACmpt.replace(0, AReCmpt);
            ACmpt.replace(1, AImCmpt);

            // Create component equations
            fvScalarMatrix AReCmptEqn
            (
                AReCmpt,
                AReEqn.dimensions()
            );
            fvScalarMatrix AImCmptEqn
            (
                AImCmpt,
                AImEqn.dimensions()
            );

            // Copy symmetric matrix coeffs and source
            AReCmptEqn.diag() = AReEqn.diag();
            AReCmptEqn.upper() = AReEqn.upper();
            AReCmptEqn.internalCoeffs() = AReEqn.internalCoeffs().component(cmpt);
            AReCmptEqn.boundaryCoeffs() = AReEqn.boundaryCoeffs().component(cmpt);
            AReCmptEqn.source() = AReEqn.source().component(cmpt);
            AImCmptEqn.diag() = AImEqn.diag();
            AImCmptEqn.upper() = AImEqn.upper();
            AImCmptEqn.internalCoeffs() = AImEqn.internalCoeffs().component(cmpt);
            AImCmptEqn.boundaryCoeffs() = AImEqn.boundaryCoeffs().component(cmpt);
            AImCmptEqn.source() = AImEqn.source().component(cmpt);

            // Add explicit gradient terms of electric scalar potential V
            AReCmptEqn -= sigma * VReGrad.component(cmpt);
            AImCmptEqn -= sigma * VImGrad.component(cmpt);

            // Prepare block system for current component
            fvBlockMatrix<vector2> AcmptEqn(ACmpt);

            // Insert equations into block matrix
            AcmptEqn.insertEquation(0, AReCmptEqn);
            AcmptEqn.insertEquation(1, AImCmptEqn);

            // Coupling coeff field
//             const scalarField alpha = omega * sigma;

            // Add coupling terms for real/imagninary ARe-AIm
            // dependence of magnetic vector potential
            AcmptEqn.insertEquationCoupling(0, 1,  omega * sigma);
            AcmptEqn.insertEquationCoupling(1, 0, -omega * sigma);

            // Solve A
            BlockSolverPerformance<vector2> AcmptEqnPerf = AcmptEqn.solve(Adict);

            // Store residuals
            scalar resAcmptRe = AcmptEqnPerf.initialResidual().component(0);
            scalar resAcmptIm = AcmptEqnPerf.initialResidual().component(1);
            resAcmpt = max(resAcmpt, max(resAcmptRe, resAcmptIm));

            // Retrieve solution
            AcmptEqn.retrieveSolution(0, AReCmpt.internalField());
            AcmptEqn.retrieveSolution(1, AImCmpt.internalField());

            // Correct component boundary conditions
            AReCmpt.correctBoundaryConditions();
            AImCmpt.correctBoundaryConditions();

            // Replace components of A
            ARe.replace(cmpt, AReCmpt);
            AIm.replace(cmpt, AImCmpt);
        }

        // Store max residual
        resA = resAcmpt;

        // Correct boundary conditions
        ARe.correctBoundaryConditions();
        AIm.correctBoundaryConditions();
    }
