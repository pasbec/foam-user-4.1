
    // Solve A component-wise in default region
    if (!(globalSettings.biotSavart && globalSettings.lowFrequency))
    {
        SM_REGIONSCOPE(DEFAULT);

        volScalarField& sigma = globalStorage.sigma()[DEFAULT];
        volScalarField& rMur = globalStorage.rMur()[DEFAULT];

        volVectorField& ARe = globalStorage.ARe()[DEFAULT];
        volVectorField& AIm = globalStorage.AIm()[DEFAULT];

// TEST
//         volScalarField& GRe = storage.GRe();
//         volScalarField& GIm = storage.GIm();

// TEST
//         // Get under-relaxation factor
//         const scalar AReUrf = mesh.solutionDict().equationRelaxationFactor(ARe.name());
//         const scalar AImUrf = mesh.solutionDict().equationRelaxationFactor(AIm.name());
//         const scalar GReUrf = mesh.solutionDict().fieldRelaxationFactor(GRe.name());
//         const scalar GImUrf = mesh.solutionDict().fieldRelaxationFactor(GIm.name());

// TEST
        {
            sigma.write();
            rMur.write();

            ARe.rename("AReOld");
            AIm.rename("AImOld");

            ARe.write();
            AIm.write();

            ARe.rename("ARe");
            AIm.rename("AIm");

            if (globalSettings.biotSavart)
            {
                globalStorage.A0Re()[DEFAULT].write();
                globalStorage.A0Im()[DEFAULT].write();
            }

            globalStorage.VReGrad()[DEFAULT].write();
            globalStorage.VImGrad()[DEFAULT].write();
        }

// TEST
        {
            // Eddy current density
            volVectorField jRe("jRe", omega0 * sigma * AIm);
            volVectorField jIm("jIm", -omega0 * sigma * ARe);

            if (control.meshIs3D())
            {
                jRe -= sigma * globalStorage.VReGrad()[DEFAULT];
                jIm -= sigma * globalStorage.VImGrad()[DEFAULT];
            }

            if (globalSettings.biotSavart)
            {
                jRe += omega0 * sigma * globalStorage.A0Im()[DEFAULT];
                jIm -= omega0 * sigma * globalStorage.A0Re()[DEFAULT];
            }

            jRe.write();
            jIm.write();

            volScalarField jReDiv = fvc::div(jRe);
            volScalarField jImDiv = fvc::div(jIm);

            Info<< Control::typeName << " | AEqn.H : "
                << "max(mag(jReDiv)) = " << max(mag(jReDiv)).value()
                << endl;
            Info<< Control::typeName << " | AEqn.H : "
                << "max(mag(jImDiv)) = " << max(mag(jImDiv)).value()
                << endl;

            jReDiv.write();
            jImDiv.write();
        }

// TEST
        {
            ARe.write();
            AIm.write();

            volScalarField AReDiv = fvc::div(ARe);
            volScalarField AImDiv = fvc::div(AIm);

            Info<< Control::typeName << " | AEqn.H : "
                << "max(mag(AReDiv)) = " << max(mag(AReDiv)).value()
                << endl;
            Info<< Control::typeName << " | AEqn.H : "
                << "max(mag(AImDiv)) = " << max(mag(AImDiv)).value()
                << endl;

            AReDiv.write();
            AImDiv.write();
        }

        // Interpolate reciprocal relative magnetic permeability
        surfaceScalarField rMurf("rMurf", fvc::interpolate(rMur));

        // Assemble base-matrices for A
        tmp<fvVectorMatrix> tAReEqn
        (
          - fvm::laplacian(rMurf,ARe)
        );

        tmp<fvVectorMatrix> tAImEqn
        (
          - fvm::laplacian(rMurf,AIm)
        );

        fvVectorMatrix& AReEqn = tAReEqn();
        fvVectorMatrix& AImEqn = tAImEqn();

        // Add source current or source vector potential if Biot-Savart is used
        if (!globalSettings.biotSavart)
        {
            if (Control::debug)
            {
                Info<< Control::typeName << " | AEqn.H : "
                    << "Add source current density to AEqn"
                    << endl;
            }

            AReEqn -= physicalConstant::mu0 * storage.j0Re();
            AImEqn -= physicalConstant::mu0 * storage.j0Im();
        }
        else
        {
            if (Control::debug)
            {
                Info<< Control::typeName << " | AEqn.H : "
                    << "Add impressed vector potenial to AEqn"
                    << endl;
            }

            AReEqn -= physicalConstant::mu0 * sigma * omega0 * globalStorage.A0Im()[DEFAULT];
            AImEqn += physicalConstant::mu0 * sigma * omega0 * globalStorage.A0Re()[DEFAULT];
        }

        // Add explicit gradient terms of electric scalar potential V
        if (control.meshIs3D() && !globalSettings.lowFrequency)
        {
            if (Control::debug)
            {
                Info<< Control::typeName << " | AEqn.H : "
                    << "Add gradient of electric scalar potential to AEqn"
                    << endl;
            }

            AReEqn += physicalConstant::mu0 * sigma * globalStorage.VReGrad()[DEFAULT];
            AImEqn += physicalConstant::mu0 * sigma * globalStorage.VImGrad()[DEFAULT];
        }

        // Create list with boundary types
        wordList ARePatchTypes(mesh.boundary().size());
        wordList AImPatchTypes(mesh.boundary().size());

        forAll(ARe.boundaryField(), patchI)
        {
            ARePatchTypes[patchI] = ARe.boundaryField()[patchI].type();
        }
        forAll(AIm.boundaryField(), patchI)
        {
            AImPatchTypes[patchI] = AIm.boundaryField()[patchI].type();
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
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
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
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                AIm.component(cmpt),
                AImPatchTypes
            );

            // Correct/Evaluate boundary values
            AReCmpt.correctBoundaryConditions();
            AImCmpt.correctBoundaryConditions();

            // Create component equations
            tmp<fvScalarMatrix> tAReCmptEqn
            (
                new fvScalarMatrix
                (
                    AReCmpt,
                    AReEqn.dimensions()
                )
            );
            tmp<fvScalarMatrix> tAImCmptEqn
            (
                new fvScalarMatrix
                (
                    AImCmpt,
                    AImEqn.dimensions()
                )
            );

            fvScalarMatrix& AReCmptEqn = tAReCmptEqn();
            fvScalarMatrix& AImCmptEqn = tAImCmptEqn();

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

// // TEST
//             AReCmptEqn.relax(AReUrf);
//             AImCmptEqn.relax(AImUrf);
//
//             AReCmptEqn += fvc::grad(GRe)().component(cmpt);
//             AImCmptEqn += fvc::grad(GIm)().component(cmpt);

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

            // Clear memory
            tAReCmptEqn.clear();
            tAImCmptEqn.clear();

            // Add coupling terms for real/imagninary ARe-AIm
            // dependence of magnetic vector potential
            if (!globalSettings.lowFrequency)
            {
                AcmptEqn.insertEquationCoupling(0, 1, - physicalConstant::mu0 * sigma * omega0);
                AcmptEqn.insertEquationCoupling(1, 0,   physicalConstant::mu0 * sigma * omega0);
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

        // Clear memory
        tAReEqn.clear();
        tAImEqn.clear();

// TEST
//         {
//             volScalarField rAARe("rAARe", 1.0/AReEqn.A());
//             volScalarField rAAIm("rAAIm", 1.0/AImEqn.A());
//
//             volVectorField HARe("HARe", AReEqn.H());
//             volVectorField HAIm("HAIm", AImEqn.H());
//
//             AReEqn -= physicalConstant::mu0 * sigma * omega0 * rAAIm * HAIm;
//             AImEqn += physicalConstant::mu0 * sigma * omega0 * rAARe * HARe;
//
//             AReEqn -= physicalConstant::mu0 * sigma * omega0 * rAAIm * fvc::grad(GIm);
//             AImEqn += physicalConstant::mu0 * sigma * omega0 * rAARe * fvc::grad(GRe);
//
//             AReEqn -= fvm::Sp(pow(physicalConstant::mu0 * sigma * omega0, 2) * rAAIm, ARe);
//             AImEqn -= fvm::Sp(pow(physicalConstant::mu0 * sigma * omega0, 2) * rAARe, AIm);
//         }

// TEST
//         solve(AReEqn);
//         solve(AImEqn);
//         solve(AReEqn == -GReGrad);
//         solve(AImEqn == -GImGrad);
//         solve(relax(AReEqn, AReUrf) == -fvc::grad(GRe));
//         solve(relax(AImEqn, AImUrf) == -fvc::grad(GIm));

// // TEST
//         for (int AEqnI = 0; AEqnI < 2; AEqnI++)
//         {
//             volVectorField* APtr = NULL;
//             fvVectorMatrix* AEqnPtr = NULL;
//             volScalarField* GPtr = NULL;
//             scalar AUrf = 1.0;
//             scalar GUrf = 1.0;
//
//             if (AEqnI == 0)
//             {
//                 APtr = &ARe;
//                 AEqnPtr = &AReEqn;
//                 GPtr = &GRe;
//                 AUrf = AReUrf;
//                 GUrf = GReUrf;
//             }
//             else
//             {
//                 APtr = &AIm;
//                 AEqnPtr = &AImEqn;
//                 GPtr = &GIm;
//                 AUrf = AImUrf;
//                 GUrf = GImUrf;
//             }
//
//             volVectorField& A = *APtr;
//             fvVectorMatrix& AEqn = *AEqnPtr;
//             volScalarField& G = *GPtr;
//
//             volScalarField rAA("rA" + A.name(), 1.0/AEqn.A());
//             surfaceScalarField rAAf("rA" + A.name() + "f", fvc::interpolate(rAA));
//
//             A.storePrevIter();
//             G.storePrevIter();
//
// // TODO: 2D case? H-operator wrongly evaluates valid components!
//             A = rAA*AEqn.H();
//
//             surfaceScalarField phiA = fvc::interpolate(A) & mesh.Sf();
//
// //             for (int i = 0; i < 3; i++)
//             {
//                 fvScalarMatrix GEqn
//                 (
//                     min(rMur) * fvm::laplacian(rAAf, G) == min(rMur) * fvc::div(phiA)
//                 );
//
//                 GEqn.solve();
//             }
//
//             G.write();
//
//             G.relax(GUrf);
//
//             A = AUrf*(A - rAA*fvc::grad(G)) + (1 - AUrf)*A.prevIter();
//             A.correctBoundaryConditions();
//         }

// TEST
//         surfaceTensorField gradARef =
//             linearInterpolate
//             (
//                 fvc::reconstruct
//                 (
//                     AReEqn.flux()/ownerInterpolate(rMur),
//                     AReEqn.jumpFlux()/neighbourInterpolate(rMur)
//                 )
//             );
//         surfaceTensorField gradAImf =
//             linearInterpolate
//             (
//                 fvc::reconstruct
//                 (
//                     AImEqn.flux()/ownerInterpolate(rMur),
//                     AImEqn.jumpFlux()/neighbourInterpolate(rMur)
//                 )
//             );
//
//         storage.gradAReNPhi() = -gradARef & mesh.Sf();
//         storage.gradAImNPhi() = -gradAImf & mesh.Sf();
//
//         storage.gradAReNPhi() -= (storage.gradAReNPhi() & mesh.Sf()) * mesh.Sf() / pow(mesh.magSf(), 2);
//         storage.gradAImNPhi() -= (storage.gradAImNPhi() & mesh.Sf()) * mesh.Sf() / pow(mesh.magSf(), 2);

// TEST
        for (int AEqnI = 0; AEqnI < 2; AEqnI++)
        {
            volVectorField* APtr = NULL;
            surfaceVectorField* gradNPhiPtr = NULL;

            if (AEqnI == 0)
            {
                APtr = &globalStorage.ARe()[DEFAULT];
                gradNPhiPtr = &storage.gradAReNPhi();
            }
            else
            {
                APtr = &globalStorage.AIm()[DEFAULT];
                gradNPhiPtr = &storage.gradAImNPhi();
            }

            volVectorField& A = *APtr;
            surfaceVectorField& gradNPhi = *gradNPhiPtr;

            surfaceTensorField gradNf =
                fvc::interpolate
                (
                    fvc::grad(linearInterpolate(A)),
                    "interpolate(grad(" + A.name() + "))"
                );

            gradNPhi = -gradNf & mesh.Sf();

            gradNPhi += tr(gradNf) * mesh.Sf();
        }
    }
