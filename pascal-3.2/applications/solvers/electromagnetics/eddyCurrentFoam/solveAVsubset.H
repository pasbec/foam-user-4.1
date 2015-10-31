
    // Solve A
    {
        // Coupling terms of electric scalar potential V
        // in equation for magnetic vector potential A
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

            // Explicit gradient terms of electric scalar potential V
            AxReEqn -= sigma * VReGrad.component(0) * indicatorVc;
            AxImEqn -= sigma * VImGrad.component(0) * indicatorVc;
            AyReEqn -= sigma * VReGrad.component(1) * indicatorVc;
            AyImEqn -= sigma * VImGrad.component(1) * indicatorVc;
            AzReEqn -= sigma * VReGrad.component(2) * indicatorVc;
            AzImEqn -= sigma * VImGrad.component(2) * indicatorVc;
        }

        // References to internal alpha field
        const scalarField& alphaIn = alpha.internalField();

        // Prepare block systems for A
        fvBlockMatrix<vector2> AxEqn(Ax);
        fvBlockMatrix<vector2> AyEqn(Ay);
        fvBlockMatrix<vector2> AzEqn(Az);

        // Insert equations for A
        AxEqn.insertEquation(0, AxReEqn);
        AxEqn.insertEquation(1, AxImEqn);
        AyEqn.insertEquation(0, AyReEqn);
        AyEqn.insertEquation(1, AyImEqn);
        AzEqn.insertEquation(0, AzReEqn);
        AzEqn.insertEquation(1, AzImEqn);

        // Coupling terms for real/imagninary ARe-AIm
        // dependence of magnetic vector potential
        {
            AxEqn.insertEquationCoupling(0, 1,  alphaIn);
            AxEqn.insertEquationCoupling(1, 0, -alphaIn);
            AyEqn.insertEquationCoupling(0, 1,  alphaIn);
            AyEqn.insertEquationCoupling(1, 0, -alphaIn);
            AzEqn.insertEquationCoupling(0, 1,  alphaIn);
            AzEqn.insertEquationCoupling(1, 0, -alphaIn);
        }

        // Solve A
        AxEqn.solve();
        AyEqn.solve();
        AzEqn.solve();

        // Retrieve solution
        AxEqn.retrieveSolution(0, AxRe.internalField());
        AxEqn.retrieveSolution(1, AxIm.internalField());
        AyEqn.retrieveSolution(0, AyRe.internalField());
        AyEqn.retrieveSolution(1, AyIm.internalField());
        AzEqn.retrieveSolution(0, AzRe.internalField());
        AzEqn.retrieveSolution(1, AzIm.internalField());

        // Replace components
        ARe.replace(0,AxRe);
        ARe.replace(1,AyRe);
        ARe.replace(2,AzRe);
        AIm.replace(0,AxIm);
        AIm.replace(1,AyIm);
        AIm.replace(2,AzIm);

        ARe.correctBoundaryConditions();
        AIm.correctBoundaryConditions();
    }

//     // Solve A
//     {
//         // Coupling terms of electric scalar potential V
//         // in equation for magnetic vector potential A
//         {
//             // Get relaxation factor
//             scalar relax = AVmethodDict.lookupOrDefault<scalar>("relax",1.0);
//
//             if (relax < (1.0-SMALL))
//             {
//                 // Relaxation of gradient of V
//                 VReGrad.relax(relax);
//                 VImGrad.relax(relax);
//             }
//
//             // Explicit gradient terms of electric scalar potential V
//             AReEqn -= sigma * VReGrad * indicatorVc;
//             AImEqn -= sigma * VImGrad * indicatorVc;
//         }
//
//         // References to internal alpha field
//         const scalarField& alphaIn = alpha.internalField();
//
//         // Prepare block system for A
//         fvBlockMatrix<vector6> AEqn(A);
//
//         // Insert equations for A
//         AEqn.insertEquation(0, AReEqn);
//         AEqn.insertEquation(3, AImEqn);
//
//         // Coupling terms for real/imagninary ARe-AIm
//         // dependence of magnetic vector potential
//         {
//             AEqn.insertEquationCoupling(0, 3,  alphaIn);
//             AEqn.insertEquationCoupling(3, 0, -alphaIn);
//             AEqn.insertEquationCoupling(1, 4,  alphaIn);
//             AEqn.insertEquationCoupling(4, 1, -alphaIn);
//             AEqn.insertEquationCoupling(2, 5,  alphaIn);
//             AEqn.insertEquationCoupling(5, 2, -alphaIn);
//         }
//
//         // Solve A
//         AEqn.solve();
//
//         // Retrieve solution
//         AEqn.retrieveSolution(0, ARe.internalField());
//         AEqn.retrieveSolution(3, AIm.internalField());
//
//         ARe.correctBoundaryConditions();
//         AIm.correctBoundaryConditions();
//
//         // TODO: Remove after debug
//         {
//             vector6Field AResIn = AEqn.residual(A.internalField());
//
//             volVectorField AReRes
//             (
//                 IOobject
//                 (
//                     "AReRes",
//                     mesh.time().timeName(),
//                     mesh,
//                     IOobject::NO_READ,
//                     IOobject::NO_WRITE
//                 ),
//                 mesh,
//                 dimensioned<vector>(word(),dimless,pTraits<vector>::zero),
//                 fixedValueFvPatchField<vector>::typeName
//             );
//             volVectorField AImRes
//             (
//                 IOobject
//                 (
//                     "AImRes",
//                     mesh.time().timeName(),
//                     mesh,
//                     IOobject::NO_READ,
//                     IOobject::NO_WRITE
//                 ),
//                 mesh,
//                 dimensioned<vector>(word(),dimless,pTraits<vector>::zero),
//                 fixedValueFvPatchField<vector>::typeName
//             );
//
//             AReRes.internalField().replace(0, AResIn.component(0));
//             AReRes.internalField().replace(1, AResIn.component(1));
//             AReRes.internalField().replace(2, AResIn.component(2));
//             AImRes.internalField().replace(0, AResIn.component(4));
//             AImRes.internalField().replace(1, AResIn.component(5));
//             AImRes.internalField().replace(2, AResIn.component(6));
//
//             AReRes.write();
//             AImRes.write();
//         }
//     }

    // Solve V
    {
        // Create interpolated surface normal components of A
        surfaceVectorField sARe("sARe", fvc::interpolate(ARe));
        surfaceVectorField sAIm("sAIm", fvc::interpolate(AIm));

        // Create local dictionory for magnetic properties
        // TODO: This is currently necessary for
        //       electricPotentialGrad
        IOdictionary cMagneticProperties
        (
            IOobject
            (
                "magneticProperties",
                cMesh.time().constant(),
                cMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            magneticProperties
        );

        // Create local magnetic vector potential A
        volVectorField cARe
        (
            IOobject
            (
                "ARe",
                cMesh.time().timeName(),
                cMesh,
                IOobject::NO_READ
            ),
            cMesh,
            dimensioned<vector>(word(), ARe.dimensions(), vector::zero),
            fixedValueFvPatchField<vector>::typeName
        );
        volVectorField cAIm
        (
            IOobject
            (
                "AIm",
                cMesh.time().timeName(),
                cMesh,
                IOobject::NO_READ
            ),
            cMesh,
            dimensioned<vector>(word(), AIm.dimensions(), vector::zero),
            fixedValueFvPatchField<vector>::typeName
        );

        // Create local electric scalar potential V
        volScalarField cVRe
        (
            IOobject
            (
                "VRe",
                cMesh.time().timeName(),
                cMesh,
                IOobject::NO_READ
            ),
            cMesh,
            dimensioned<scalar>(word(), VRe.dimensions(), 0.0),
            fixedGradientFvPatchField<scalar>::typeName
        );
        volScalarField cVIm
        (
            IOobject
            (
                "VIm",
                cMesh.time().timeName(),
                cMesh,
                IOobject::NO_READ
            ),
            cMesh,
            dimensioned<scalar>(word(), VIm.dimensions(), 0.0),
            fixedGradientFvPatchField<scalar>::typeName
        );

        // Create local conductivity
        volScalarField cSigma
        (
            IOobject
            (
                "sigma",
                cMesh.time().timeName(),
                cMesh,
                IOobject::NO_READ
            ),
            cMesh,
            dimensionedScalar(word(), sigma.dimensions(), 0.0),
            zeroGradientFvPatchField<scalar>::typeName
        );

        // Internal fields
        forAll (cMesh.C(), cellI)
        {
            // Map internal values for A
            cARe.internalField()[cellI] = ARe.internalField()[cMeshCellMap[cellI]];
            cAIm.internalField()[cellI] = AIm.internalField()[cMeshCellMap[cellI]];

            // Map internal values for V
            cVRe.internalField()[cellI] = VRe.internalField()[cMeshCellMap[cellI]];
            cVIm.internalField()[cellI] = VIm.internalField()[cMeshCellMap[cellI]];

            // Map internal values for sigma
            cSigma.internalField()[cellI] = sigma.internalField()[cMeshCellMap[cellI]];
        }

        // Boundary fields
        forAll(cMesh.boundaryMesh(), patchI)
        {
            polyPatch patch = cMesh.boundaryMesh()[patchI];

            // Old internal faces
            if (patch.name() == "oldInternalFaces")
            {
                // Boundary values for A
                {
                    vectorField& cAReP = cARe.boundaryField()[patchI];
                    vectorField& cAImP = cAIm.boundaryField()[patchI];

                    forAll(patch, facei)
                    {
                        label faceI = patch.start() + facei;

                        cAReP[facei] = sARe[cMeshFaceMap[faceI]];
                        cAImP[facei] = sAIm[cMeshFaceMap[faceI]];
                    }
                }

                // Boundary values for V
                {
                    if
                    (
                        cVRe.boundaryField()[patchI].type()
                    == fixedGradientFvPatchField<scalar>::typeName
                    )
                    {
                        fixedGradientFvPatchField<scalar>& cVReP =
                            refCast<fixedGradientFvPatchField<scalar> >
                            (
                                cVRe.boundaryField()[patchI]
                            );

                        // Set correct gradient
                        cVReP.gradient() =  omega.value()
                        * (cAIm.boundaryField()[patchI] & cVReP.patch().nf());

                        // Evaluate values
                        cVReP.evaluate();
                    }

                    if
                    (
                        cVIm.boundaryField()[patchI].type()
                    == fixedGradientFvPatchField<scalar>::typeName
                    )
                    {
                        fixedGradientFvPatchField<scalar>& cVImP =
                            refCast<fixedGradientFvPatchField<scalar> >
                            (
                                cVIm.boundaryField()[patchI]
                            );

                        // Set correct gradient
                        cVImP.gradient() = -omega.value()
                        * (cARe.boundaryField()[patchI] & cVImP.patch().nf());

                        // Evaluate values
                        cVImP.evaluate();
                    }
                }

                // Boundary values for sigma
                {
                    fvPatchField<scalar>& cSigmaP = cSigma.boundaryField()[patchI];

                    // Evaluate values
                    cSigmaP.evaluate();
                }
            }

            // TODO: What about other patches?
        }

        // Correct boundary conditions for A
        cARe.correctBoundaryConditions();
        cAIm.correctBoundaryConditions();

        // Correct boundary conditions for V
        cVRe.correctBoundaryConditions();
        cVIm.correctBoundaryConditions();

        // Correct boundary conditions for sigma
        cSigma.correctBoundaryConditions();

        // Prepare V-equations
        fvScalarMatrix cVReEqn(fvm::laplacian(cSigma, cVRe, "laplacian(sigma,VRe)"));
        fvScalarMatrix cVImEqn(fvm::laplacian(cSigma, cVIm, "laplacian(sigma,VIm)"));

        volVectorField cSigmaGrad = fvc::grad(cSigma);

        // Explicit source terms due to variable sigma
        // in conducting region with implied condition div(A)=0
        volScalarField cVReEqnSource =   omega * (cSigmaGrad & cAIm);
        volScalarField cVImEqnSource = - omega * (cSigmaGrad & cARe);

//         // TEST
//         cVReEqnSource +=   omega * cSigma * fvc::div(cAIm);
//         cVImEqnSource += - omega * cSigma * fvc::div(cARe);

        // Check/Correct Poisson compatibility
        // TEST
        {
            // Source sums
            scalar cVReEqnSumSource = 0.0;
            scalar cVImEqnSumSource = 0.0;

            // Flux sums
            scalar cVReEqnSumFlux = 0.0;
            scalar cVImEqnSumFlux = 0.0;

            // Correct source integrals
            // TEST
            // FIXME: Parallel?
            {
                // Global cell size
                label nCellsGlobal = returnReduce(cMesh.nCells(), sumOp<label>());

                // Global source integrals
                cVReEqnSumSource = gSum(cVReEqnSource.internalField() * cMesh.V());
                cVImEqnSumSource = gSum(cVImEqnSource.internalField() * cMesh.V());

                // Force global source integral to zero
                cVReEqnSource.internalField() -= cVReEqnSumSource / nCellsGlobal / cMesh.V();
                cVImEqnSource.internalField() -= cVImEqnSumSource / nCellsGlobal / cMesh.V();
            }

            // Source integrals
            cVReEqnSumSource = gSum(cVReEqnSource.internalField() * cMesh.V());
            cVImEqnSumSource = gSum(cVImEqnSource.internalField() * cMesh.V());

            // Boundary integrals
            {
                forAll(cMesh.boundaryMesh(), patchI)
                {
                    const polyPatch& patch = cMesh.boundaryMesh()[patchI];

                    const vectorField& SfP = cMesh.Sf().boundaryField()[patchI];

                    const vectorField& cAReP = cARe.boundaryField()[patchI];
                    const vectorField& cAImP = cAIm.boundaryField()[patchI];

                    // Old internal faces
                    if (patch.name() == "oldInternalFaces")
                    {
                        cVReEqnSumFlux += sum(omega.value() * (cAReP & SfP));
                        cVImEqnSumFlux += sum(omega.value() * (cAImP & SfP));
                    }

                    // TODO: What about other patches?
                }

                reduce(cVReEqnSumFlux, sumOp<scalar>());
                reduce(cVImEqnSumFlux, sumOp<scalar>());
            }

            Info << "VRe Poisson equation compatibility check (Source, Flux):"
                << "(" << cVReEqnSumSource << ", " << cVReEqnSumFlux << ")" << endl;
            Info << "VIm Poisson equation compatibility check (Source, Flux):"
                << "(" << cVImEqnSumSource << ", " << cVImEqnSumFlux << ")" << endl;
        }

        // Add source terms
        cVReEqn -= cVReEqnSource;
        cVImEqn -= cVImEqnSource;

//         // Reference point
//         // FIXME: Refpoint problems!
//         {
//             // Reference values in conducting region
//             scalar cRefValueVRe = 0.0;
//             scalar cRefValueVIm = 0.0;
//
//             // Reference cells in conducting region
//             cVReEqn.setReference(csRefCell, cRefValueVRe);
//             cVImEqn.setReference(csRefCell, cRefValueVIm);
//         }

//         // Pseudo time-stepping based on Euler-explicit
//         // TEST
//         {
//             // Grab parts of V-equations
//             scalarField& cVReEqnD = cVReEqn.diag();
//             scalarField& cVReEqnS = cVReEqn.source();
//             scalarField& cVImEqnD = cVImEqn.diag();
//             scalarField& cVImEqnS = cVImEqn.source();
//
// //             // Fictive time step
// //             // NOTE: Minus sign as laplacian is in matrix with plus sign
// //             scalarField V = cMesh.V();
// //             scalar pseudoDeltaT = 1.0;
// //             cVReEqnD -= V/pseudoDeltaT;
// //             cVImEqnD -= V/pseudoDeltaT;
// //             cVReEqnS -= V/pseudoDeltaT * cVRe.internalField();
// //             cVImEqnS -= V/pseudoDeltaT * cVIm.internalField();
//
//             // Arbitrary damp factor based on diagonal value
//             scalar damp = 1.0e-04;
//             cVReEqnD += damp * cVReEqnD;
//             cVImEqnD += damp * cVImEqnD;
//             cVReEqnS += damp * cVReEqnD * cVRe.internalField();
//             cVImEqnS += damp * cVImEqnD * cVIm.internalField();
//         }

        // NOTE NOTE NOTE: Currently deflated solver is necessary as matrix is kept singluar!

        // Restore solution level
        cVRe += VReAverage;
        cVIm += VImAverage;

        // Solve V-equations
        cVReEqn.solve();
        cVImEqn.solve();

        // Correct boundary conditions
        cVRe.correctBoundaryConditions();
        cVIm.correctBoundaryConditions();

        // Calculate average
        VReAverage = average(cVRe);
        VImAverage = average(cVIm);

        // Normalize solution level
        cVRe -= VReAverage;
        cVIm -= VImAverage;

//         // Dirichlet correction
//         // TEST
//         {
//             // Create copy of electric scalar potential V
//             // based on fixed value boundary conditions
//             volScalarField cVReCorr
//             (
//                 IOobject
//                 (
//                     "VReCorr",
//                     cMesh.time().timeName(),
//                     cMesh,
//                     IOobject::NO_READ
//                 ),
//                 cMesh,
//                 dimensioned<scalar>(word(), VRe.dimensions(), 0.0),
//                 fixedValueFvPatchField<scalar>::typeName
//             );
//             volScalarField cVImCorr
//             (
//                 IOobject
//                 (
//                     "VImCorr",
//                     cMesh.time().timeName(),
//                     cMesh,
//                     IOobject::NO_READ
//                 ),
//                 cMesh,
//                 dimensioned<scalar>(word(), VIm.dimensions(), 0.0),
//                 fixedValueFvPatchField<scalar>::typeName
//             );
//
//             // Copy Neumann-Solution including evaluated boundary values
//             cVReCorr == cVRe;
//             cVImCorr == cVIm;
//
//             // Prepare V-equations based on fixed value boundary conditions
//             fvScalarMatrix cVReCorrEqn
//             (
//                 fvm::laplacian(cSigma, cVReCorr, "laplacian(sigma,VRe)")
//             );
//             fvScalarMatrix cVImCorrEqn
//             (
//                 fvm::laplacian(cSigma, cVImCorr, "laplacian(sigma,VIm)")
//             );
//
//             // Add source to matrix
//             cVReCorrEqn -= cVReEqnSource;
//             cVImCorrEqn -= cVImEqnSource;
//
//             // Solve V-equations based on fixed value boundary conditions
//             cVReCorrEqn.solve(cMesh.solutionDict().solver("VReCorr"));
//             cVImCorrEqn.solve(cMesh.solutionDict().solver("VImCorr"));
//
//             // Correct Neumann-Solution
//             cVRe == cVReCorr;
//             cVIm == cVImCorr;
//         }

        // Calculate gradient
        volVectorField cVReGrad = fvc::grad(cVRe);
        volVectorField cVImGrad = fvc::grad(cVIm);

        // Reverse map V and its gradient
        forAll (cMesh.C(), cellI)
        {
            VRe[cMeshCellMap[cellI]] = cVRe[cellI];
            VIm[cMeshCellMap[cellI]] = cVIm[cellI];

            VReGrad[cMeshCellMap[cellI]] = cVReGrad[cellI];
            VImGrad[cMeshCellMap[cellI]] = cVImGrad[cellI];
        }
    }
