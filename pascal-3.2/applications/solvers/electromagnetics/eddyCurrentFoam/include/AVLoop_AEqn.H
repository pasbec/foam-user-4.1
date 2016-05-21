
    // Solve A component-wise in default region
    {
        SM_REGIONSCOPE(DEFAULT);

        volScalarField& sigma = globalStorage.sigma()[DEFAULT];

        volVectorField& ARe = globalStorage.ARe()[DEFAULT];
        volVectorField& AIm = globalStorage.AIm()[DEFAULT];

        volVectorField& VReGrad = globalStorage.VReGrad()[DEFAULT];
        volVectorField& VImGrad = globalStorage.VImGrad()[DEFAULT];

        // Assemble base-matrices for A
        fvVectorMatrix AReEqn
        (
            physicalConstant::rMu0 * fvm::laplacian(ARe)
//        - physicalConstant::rMu0 * fvc::grad(fvc::div(ARe))
//        + sigma * omega0 * AIm
//        - sigma * VReGrad
//        + j0Re
//        + sigma * omega0 * A0Im
        );
        fvVectorMatrix AImEqn
        (
            physicalConstant::rMu0 * fvm::laplacian(AIm)
//        - physicalConstant::rMu0 * fvc::grad(fvc::div(AIm))
//        - sigma * omega0 * ARe
//        - sigma * VImGrad
//        + j0Im
//        - sigma * omega0 * A0Re
        );

        // Add source current or source vector potential if Biot-Savart is used
        if (!globalSettings.biotSavart)
        {
            AReEqn += storage.j0Re();
            AImEqn += storage.j0Im();
        }
        else
        {
            AReEqn += sigma * omega0 * globalStorage.A0Im()[DEFAULT];
            AImEqn -= sigma * omega0 * globalStorage.A0Re()[DEFAULT];
        }

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

        // Loop over components
        for (direction cmpt = 0; cmpt < vector::nComponents; cmpt++)
        {
            // Only solve for valid directions
            if(control.solutionDir()[cmpt] < 1) continue;

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

            // Add coupling terms for real/imagninary ARe-AIm
            // dependence of magnetic vector potential
            AcmptEqn.insertEquationCoupling(0, 1,  omega0 * sigma);
            AcmptEqn.insertEquationCoupling(1, 0, -omega0 * sigma);

            // Solve A
            // and prevent segfault for "0"-Solution
            scalar preResidual= gMax(mag(AcmptEqn.residual(ACmpt, AcmptEqn.source())));
            if (preResidual > 0)
            {
                AcmptEqn.solve(control.subDict(Region::DEFAULT, ACmpt.name()));
            }

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

        // Correct boundary conditions
        ARe.correctBoundaryConditions();
        AIm.correctBoundaryConditions();
    }
