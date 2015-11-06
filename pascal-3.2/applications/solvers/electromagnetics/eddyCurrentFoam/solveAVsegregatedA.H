
    // Solve A component-wise
    {
        // Relax gradient of V
        {
            // Store old gradient
            VReGrad.storePrevIter();
            VImGrad.storePrevIter();

            // Get relaxation factor
            scalar relax = AVmethodDict.lookupOrDefault<scalar>("relax",1.0);

            if (relax < (1.0-SMALL))
            {
                // Relaxation of gradient of V
                VReGrad.relax(relax);
                VImGrad.relax(relax);
            }
        }

        // Loop over components
        forAll(vector::one, cI)
        {
            // Create fields
            volScalarField AReCmpt
            (
                IOobject
                (
                    "A" + word(vector::componentNames[cI]) + "Re",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                ARe.component(cI),
                ARePatchTypes
            );
            volScalarField AImCmpt
            (
                IOobject
                (
                    "A" + word(vector::componentNames[cI]) + "Im",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                AIm.component(cI),
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
                    "A" + word(vector::componentNames[cI]),
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
            AReCmptEqn.internalCoeffs() = AReEqn.internalCoeffs().component(cI);
            AReCmptEqn.boundaryCoeffs() = AReEqn.boundaryCoeffs().component(cI);
            AReCmptEqn.source() = AReEqn.source().component(cI);
            AImCmptEqn.diag() = AImEqn.diag();
            AImCmptEqn.upper() = AImEqn.upper();
            AImCmptEqn.internalCoeffs() = AImEqn.internalCoeffs().component(cI);
            AImCmptEqn.boundaryCoeffs() = AImEqn.boundaryCoeffs().component(cI);
            AImCmptEqn.source() = AImEqn.source().component(cI);

            // Add explicit gradient terms of electric scalar potential V
            AReCmptEqn -= sigma * VReGrad.component(cI) * indicatorVc;
            AImCmptEqn -= sigma * VImGrad.component(cI) * indicatorVc;

            // References to internal alpha field
            const scalarField& alphaIn = alpha.internalField();

            // Prepare block system for current component
            fvBlockMatrix<vector2> AcmptEqn(ACmpt);

            // Insert equations into block matrix
            AcmptEqn.insertEquation(0, AReCmptEqn);
            AcmptEqn.insertEquation(1, AImCmptEqn);

            // Add coupling terms for real/imagninary ARe-AIm
            // dependence of magnetic vector potential
            AcmptEqn.insertEquationCoupling(0, 1,  alphaIn);
            AcmptEqn.insertEquationCoupling(1, 0, -alphaIn);

            // Solve A
            AcmptEqn.solve();

            // Retrieve solution
            AcmptEqn.retrieveSolution(0, AReCmpt.internalField());
            AcmptEqn.retrieveSolution(1, AImCmpt.internalField());

            // Correct component boundary conditions
            AReCmpt.correctBoundaryConditions();
            AImCmpt.correctBoundaryConditions();

            // Replace components of A
            ARe.replace(cI,AReCmpt);
            AIm.replace(cI,AImCmpt);
        }

        // Correct boundary conditions
        ARe.correctBoundaryConditions();
        AIm.correctBoundaryConditions();
    }
