
    // Solve A component-wise in default region
    if (!(globalSettings.biotSavart && globalSettings.lowFrequency))
    {
        SM_REGIONSCOPE(DEFAULT);

        volScalarField& sigma = globalStorage.sigma()[DEFAULT];
        volScalarField& mur = globalStorage.mur()[DEFAULT];
        volScalarField mu = physicalConstant::mu0 * mur;

        volVectorField& ARe = globalStorage.ARe()[DEFAULT];
        volVectorField& AIm = globalStorage.AIm()[DEFAULT];

// TEST
        volScalarField rMu("rMu", 1.0/mu);
        surfaceVectorField gradAReNPhi("gradAReNPhi", fvc::interpolate(fvc::grad(ARe)) & mesh.Sf());
        surfaceVectorField gradAImNPhi("gradAImNPhi", fvc::interpolate(fvc::grad(AIm)) & mesh.Sf());

        // Assemble base-matrices for A
        fvVectorMatrix AReEqn
        (
            fvm::laplacian(ARe)
//             physicalConstant::rMu0/mur * fvm::laplacian(ARe)
//        - physicalConstant::rMu0/mur * fvc::grad(fvc::div(ARe))
//        + (fvc::curl(ARe) ^ fvc::grad(volScalarField("rMu", physicalConstant::rMu0/mur)))
//        + sigma * omega0 * AIm
//        - sigma * VReGrad()[DEFAULT]
//        + j0Re
//        + sigma * omega0 * A0Im
        );
        fvVectorMatrix AImEqn
        (
            fvm::laplacian(AIm)
//             physicalConstant::rMu0/mur * fvm::laplacian(AIm)
//        - physicalConstant::rMu0/mur * fvc::grad(fvc::div(AIm))
//        + (fvc::curl(AIm) ^ fvc::grad(volScalarField("rMu", physicalConstant::rMu0/mur)))
//        - sigma * omega0 * ARe
//        - sigma * VImGrad()[DEFAULT]
//        + j0Im
//        - sigma * omega0 * A0Re
        );

// TEST
//         {
//             ARe.write();
//             AIm.write();
//
//             volScalarField AReDiv = fvc::div(ARe);
//             volScalarField AImDiv = fvc::div(AIm);
//
//             AReDiv.write();
//             AImDiv.write();
//
// //             AReEqn -= (1.0-mur) * fvc::grad(AReDiv);
// //             AImEqn -= (1.0-mur) * fvc::grad(AImDiv);
//
//             volVectorField muGradrMu = mu * fvc::grad(volScalarField("rMu", physicalConstant::rMu0/mur));
//
//             AReEqn += fvc::curl(ARe) ^ muGradrMu;
//             AImEqn += fvc::curl(AIm) ^ muGradrMu;
//         }

        // Add source current or source vector potential if Biot-Savart is used
        if (!globalSettings.biotSavart)
        {
            AReEqn += mu * storage.j0Re();
            AImEqn += mu * storage.j0Im();

            if (Control::debug)
            {
                Info<< Control::typeName << " | AEqn.H : "
                    << "Add source current density to AEqn"
                    << endl;
            }
        }
        else
        {
            AReEqn += mu * sigma * omega0 * globalStorage.A0Im()[DEFAULT];
            AImEqn -= mu * sigma * omega0 * globalStorage.A0Re()[DEFAULT];

            if (Control::debug)
            {
                Info<< Control::typeName << " | AEqn.H : "
                    << "Add impressed vector potenial to AEqn"
                    << endl;
            }
        }

        // Add explicit gradient terms of electric scalar potential V
        if (control.meshIs3D() && !globalSettings.lowFrequency)
        {
            AReEqn -= mu * sigma * globalStorage.VReGrad()[DEFAULT];
            AImEqn -= mu * sigma * globalStorage.VImGrad()[DEFAULT];

            if (Control::debug)
            {
                Info<< Control::typeName << " | AEqn.H : "
                    << "Add gradient of electric scalar potential to AEqn"
                    << endl;
            }
        }

// TODO Coupled and tangential magnetic bc?
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
                dimensionedVector2(word(), dimless, vector2::zero)
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

            // Prepare block system for current component
            fvBlockMatrix<vector2> AcmptEqn(ACmpt);

            // Lookup solver type for AcmptEqn
            word solverType =
                mesh.solutionDict().subDict("solvers").
                subDict(AcmptEqn.psi().name()).lookup("solver");

            // Force diag/upper/lower square for AMG (matrix for A is asymmetric)
            if (solverType == "AMG")
            {
                Field<tensor2>& AEqnD = AcmptEqn.diag().asSquare(); AEqnD = tensor2::zero;
                Field<tensor2>& AEqnU = AcmptEqn.upper().asSquare(); AEqnU = tensor2::zero;
                Field<tensor2>& AEqnL = AcmptEqn.lower().asSquare(); AEqnL = tensor2::zero;
            }

            // Insert equations into block matrix
            AcmptEqn.insertEquation(0, AReCmptEqn);
            AcmptEqn.insertEquation(1, AImCmptEqn);

            // Add coupling terms for real/imagninary ARe-AIm
            // dependence of magnetic vector potential
            if (!globalSettings.lowFrequency)
            {
                AcmptEqn.insertEquationCoupling(0, 1,  mu * sigma * omega0);
                AcmptEqn.insertEquationCoupling(1, 0, -mu * sigma * omega0);
            }

            // Solve A but prevent segfault for homogenous solution
            scalar preResidual = mag(gSum(cmptMag(AcmptEqn.residual(ACmpt, AcmptEqn.source()))));
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
