
Info << "DEBUG: Preparations" << endl;

    // Calculate msigma
    volScalarField msigma
    (
        "msigma",
        sigma
    );

    // Manipulate msigma:
    //
    //  We will replace "zero"-values with very small
    //  values at all faces in order to avoid a singular
    //  solution
    {
        // Parameters
        scalar msigmaMin = 1.0;

        // Manipulate internal faces
        forAll(msigma.internalField(), cellI)
        {
            if
            (msigma.internalField()[cellI] < msigmaMin)
            {
                msigma.internalField()[cellI] = msigmaMin;
            }
        }

        // Manipulate boundary faces
        msigma.boundaryField() == msigmaMin;
    }

    // Create field for magnetic diffusivity
    volScalarField mgamma
    (
        "mgamma",
        rMu0/msigma
    );

    // Interpolate sigma to face centers
    surfaceScalarField msigmaf
    (
        "msigmaf",
        fvc::interpolate(msigma, "interpolate(sigma)")
    );

    // Create surface field for magnetic diffusivity
    surfaceScalarField mgammaf
    (
        "mgammaf",
        rMu0/msigmaf
    );

Info << "DEBUG: Calculate Ts" << endl;

    // Calculate Ts from source current
    {
        fvVectorMatrix TsReEqn
        (
            fvm::laplacian(TsRe)
          + fvc::curl(jsRe)
        );

        fvVectorMatrix TsImEqn
        (
            fvm::laplacian(TsIm)
          + fvc::curl(jsIm)
        );

        TsReEqn.solve();
        TsImEqn.solve();

        TsRe.correctBoundaryConditions();
        TsIm.correctBoundaryConditions();

        TsRe.write();
        TsIm.write();
    }

Info << "DEBUG: Define c/n regions" << endl;

    // Create set for conducting cells
    labelHashSet cCellsHashSet;
    forAll(mesh.cells(), cellI)
    {
        if (sigma.internalField()[cellI] > SMALL)
        {
            cCellsHashSet.insert(cellI);
        }
    }
    labelList cCells = cCellsHashSet.toc();

    // Create set for non-conducting cells
    labelHashSet nCellsHashSet;
    forAll(mesh.cells(), cellI)
    {
        if (sigma.internalField()[cellI] < SMALL)
        {
            nCellsHashSet.insert(cellI);
        }
    }
    labelList nCells = nCellsHashSet.toc();

Info << "DEBUG: Start in conducting region" << endl;

    {
        // Create submesh
        fvMeshSubset cMeshSubSet
        (
            IOobject
            (
                "conductor",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ
            ),
            mesh
        );

        // Initialize large cell subset
        cMeshSubSet.setLargeCellSubset(cCellsHashSet);

        // Create sub-mesh
        fvMesh& cMesh = cMeshSubSet.subMesh();

        // References for maps
        const labelList& cCellMap = cMeshSubSet.cellMap();

Info << "DEBUG (c): Field creation" << endl;

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
            dimensionedScalar(word(), sigma.dimensions(), pTraits<scalar>::zero),
            zeroGradientFvPatchField<scalar>::typeName
        );

//         volScalarField cSRe
//         (
//             IOobject
//             (
//                 "SRe",
//                 cMesh.time().timeName(),
//                 cMesh,
//                 IOobject::NO_READ
//             ),
//             cMesh,
//             dimensionedScalar(word(), SRe.dimensions(), pTraits<scalar>::zero),
//             fixedGradientFvPatchField<scalar>::typeName
//         );
//
//         volScalarField cSIm
//         (
//             IOobject
//             (
//                 "SIm",
//                 cMesh.time().timeName(),
//                 cMesh,
//                 IOobject::NO_READ
//             ),
//             cMesh,
//             dimensionedScalar(word(), SIm.dimensions(), pTraits<scalar>::zero),
//             fixedGradientFvPatchField<scalar>::typeName
//         );

        volScalarField cSRe
        (
            IOobject
            (
                "SRe",
                cMesh.time().timeName(),
                cMesh,
                IOobject::NO_READ
            ),
            cMesh,
            dimensionedScalar(word(), SRe.dimensions(), pTraits<scalar>::zero),
            zeroGradientFvPatchField<scalar>::typeName
        );

        volScalarField cSIm
        (
            IOobject
            (
                "SIm",
                cMesh.time().timeName(),
                cMesh,
                IOobject::NO_READ
            ),
            cMesh,
            dimensionedScalar(word(), SIm.dimensions(), pTraits<scalar>::zero),
            zeroGradientFvPatchField<scalar>::typeName
        );

        volVectorField cTRe
        (
            IOobject
            (
                "TRe",
                cMesh.time().timeName(),
                cMesh,
                IOobject::NO_READ
            ),
            cMesh,
            dimensioned<vector>(word(), TRe.dimensions(), pTraits<vector>::zero),
            newDirectionMixedFvPatchField<vector>::typeName
        );

        volVectorField cTIm
        (
            IOobject
            (
                "TIm",
                cMesh.time().timeName(),
                cMesh,
                IOobject::NO_READ
            ),
            cMesh,
            dimensionedVector(word(), TIm.dimensions(), vector::zero),
            newDirectionMixedFvPatchField<vector>::typeName
        );

        volVector6Field cT
        (
            IOobject
            (
                "T",
                cMesh.time().timeName(),
                cMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            cMesh,
            dimensionedVector6(word(),dimless,vector6::zero)
        );

        volVectorField cTsRe
        (
            IOobject
            (
                "TsRe",
                cMesh.time().timeName(),
                cMesh,
                IOobject::NO_READ
            ),
            cMesh,
            dimensioned<vector>(word(), TsRe.dimensions(), pTraits<vector>::zero),
            zeroGradientFvPatchField<vector>::typeName
        );

        volVectorField cTsIm
        (
            IOobject
            (
                "TsIm",
                cMesh.time().timeName(),
                cMesh,
                IOobject::NO_READ
            ),
            cMesh,
            dimensionedVector(word(), TsIm.dimensions(), pTraits<vector>::zero),
            zeroGradientFvPatchField<vector>::typeName
        );

        volVectorField cjRe
        (
            IOobject
            (
                "jRe",
                cMesh.time().timeName(),
                cMesh,
                IOobject::NO_READ
            ),
            cMesh,
            dimensioned<vector>(word(), jRe.dimensions(), pTraits<vector>::zero),
            calculatedFvPatchField<vector>::typeName
        );

        volVectorField cjIm
        (
            IOobject
            (
                "jIm",
                cMesh.time().timeName(),
                cMesh,
                IOobject::NO_READ
            ),
            cMesh,
            dimensionedVector(word(), jIm.dimensions(), pTraits<vector>::zero),
            calculatedFvPatchField<vector>::typeName
        );

Info << "DEBUG (c): Boundary conditions" << endl;

        // Boundary conditions
        forAll(cMesh.boundary(), patchI)
        {
//             if
//             (
//                 cSRe.boundaryField()[patchI].type()
//              == fixedGradientFvPatchField<scalar>::typeName
//             )
//             {
//                 fixedGradientFvPatchField<scalar>& cSReP =
//                     refCast<fixedGradientFvPatchField<scalar> >
//                     (
//                         cSRe.boundaryField()[patchI]
//                     );
//
//                 // Set gradient
//                 cSReP.gradient() =
//                     cTRe.boundaryField()[patchI]
//                   & cMesh.boundary()[patchI].nf();
//                   // TODO: Missing gradient from non-conducting region
//             }
//
//             if
//             (
//                 cSIm.boundaryField()[patchI].type()
//              == fixedGradientFvPatchField<scalar>::typeName
//             )
//             {
//                 fixedGradientFvPatchField<scalar>& cSImP =
//                     refCast<fixedGradientFvPatchField<scalar> >
//                     (
//                         cSIm.boundaryField()[patchI]
//                     );
//
//                 // Set gradient
//                 cSImP.gradient() =
//                     cTIm.boundaryField()[patchI]
//                   & cMesh.boundary()[patchI].nf();
//                   // TODO: Missing gradient from non-conducting region
//             }

            if
            (
                cTRe.boundaryField()[patchI].type()
             == newDirectionMixedFvPatchField<vector>::typeName
            )
            {
                newDirectionMixedFvPatchField<vector>& cTReP =
                    refCast<newDirectionMixedFvPatchField<vector> >
                    (
                        cTRe.boundaryField()[patchI]
                    );

                // Set ref value/gradient to zero
                cTReP.refValue() = vector::zero;
                cTReP.refGrad() = vector::zero;

                // Set normal to patch normal
                cTReP.nHat() = cMesh.boundary()[patchI].nf();

                // Set tangential value fraction to 1
                // to get zero tangential values
                // in order to get TRe^n=0
                cTReP.tangentialValueFraction() = 1.0;

                // Set normal value fraction to 0
                // to get zero normal derivative
                // which gives div(TRe)=0 with TRe^n=0
                cTReP.normalValueFraction() = 0.0;
            }

            if
            (
                cTIm.boundaryField()[patchI].type()
             == newDirectionMixedFvPatchField<vector>::typeName
            )
            {
                newDirectionMixedFvPatchField<vector>& cTImP =
                    refCast<newDirectionMixedFvPatchField<vector> >
                    (
                        cTIm.boundaryField()[patchI]
                    );

                // Set ref value/gradient to zero
                cTImP.refValue() = vector::zero;
                cTImP.refGrad() = vector::zero;

                // Set normal to patch normal
                cTImP.nHat() = cMesh.boundary()[patchI].nf();

                // Set tangential value fraction to 1
                // to get zero tangential values
                // in order to get TIm^n=0
                cTImP.tangentialValueFraction() = 1.0;

                // Set normal value fraction to 0
                // to get zero normal derivative
                // which gives div(TIm)=0 with TIm^n=0
                cTImP.normalValueFraction() = 0.0;
            }
        }

Info << "DEBUG (c): Mapping" << endl;

        // Map volume fields internal
        forAll (cMesh.C(), cellI)
        {
            cSigma.internalField()[cellI] = sigma.internalField()[cCellMap[cellI]];

            cSRe.internalField()[cellI] = SRe.internalField()[cCellMap[cellI]];
            cSIm.internalField()[cellI] = SIm.internalField()[cCellMap[cellI]];

            cTRe.internalField()[cellI] = TRe.internalField()[cCellMap[cellI]];
            cTIm.internalField()[cellI] = TIm.internalField()[cCellMap[cellI]];

            cTsRe.internalField()[cellI] = TsRe.internalField()[cCellMap[cellI]];
            cTsIm.internalField()[cellI] = TsIm.internalField()[cCellMap[cellI]];
        }

        // Correct boundary conditions
        cSigma.correctBoundaryConditions();

        cSRe.correctBoundaryConditions();
        cSIm.correctBoundaryConditions();

        cTRe.correctBoundaryConditions();
        cTIm.correctBoundaryConditions();

        cTsRe.correctBoundaryConditions();
        cTsIm.correctBoundaryConditions();

        // Update calculated boundary conditions
//         forAll(cMesh.boundaryMesh(), patchI)
//         {
// //             polyPatch patch = cMesh.boundaryMesh()[patchI];
// //
// //             forAll(patch, facei)
// //             {
// //                 label faceI = patch.start() + facei;
// //             }
//
//             // Use internalField for boundary field for sigma
//             cSigma.boundaryField()[patchI] = cSigma.boundaryField()[patchI].patchInternalField();
//         }

// Info << "DEBUG (c): Solve S" << endl;
//
//         volScalarField cTReDiv("TReDiv", fvc::div(cTRe)); cTReDiv.write();
//         volScalarField cTImDiv("TImDiv", fvc::div(cTIm)); cTImDiv.write();
//
//         volScalarField cTsReDiv("TsReDiv", fvc::div(cTsRe)); cTsReDiv.write();
//         volScalarField cTsImDiv("TsImDiv", fvc::div(cTsIm)); cTsImDiv.write();
//
//         fvScalarMatrix cSReEqn
//         (
//             fvm::laplacian(cSRe)
//           - fvc::div(cTRe)
//           - fvc::div(cTsRe)
//         );
//
//         fvScalarMatrix cSImEqn
//         (
//             fvm::laplacian(cSIm)
//           - fvc::div(cTIm)
//           - fvc::div(cTsIm)
//         );
//
//         // Set references
//         label refCell = mesh.findCell(point(0.00, 0.00, 0.01));
//         scalar refValue = 0.0;
//
//         cSReEqn.setReference(refCell, refValue);
//         cSImEqn.setReference(refCell, refValue);
//
//         cSReEqn.solve();
//         cSImEqn.solve();
//
//         cSRe.correctBoundaryConditions();
//         cSIm.correctBoundaryConditions();

Info << "DEBUG (c): Prepare matrices for T" << endl;

        volScalarField crMu0Sigma("rMu0Sigma", rMu0/cSigma); crMu0Sigma.write();
        volVectorField crMu0SigmaGrad("rMu0SigmaGrad", fvc::grad(rMu0/cSigma)); crMu0SigmaGrad.write();

        volVectorField cTReCurl("TReCurl", fvc::curl(cTRe)); cTReCurl.write();
        volVectorField cTImCurl("TImCurl", fvc::curl(cTIm)); cTImCurl.write();

        volVectorField cSReGrad("SReGrad", fvc::grad(cSRe)); cSReGrad.write();
        volVectorField cSImGrad("SImGrad", fvc::grad(cSIm)); cSImGrad.write();

        fvVectorMatrix cTReEqn
        (
            rMu0/cSigma * fvm::laplacian(cTRe)
          - (fvc::grad(rMu0/cSigma) ^ fvc::curl(cTRe))
//           + omega * cTIm
          + omega * cTsIm
          - omega * fvc::grad(cSIm)
        );

        fvVectorMatrix cTImEqn
        (
            rMu0/cSigma * fvm::laplacian(cTIm)
          - (fvc::grad(rMu0/cSigma) ^ fvc::curl(cTIm))
//           - omega * cTRe
          - omega * cTsRe
          + omega * fvc::grad(cSRe)
        );

Info << "DEBUG (c): Assembly block-matrix for T" << endl;

        // Prepare block system
        fvBlockMatrix<vector6> cTEqn(cT);

        // Insert equations for T
        cTEqn.insertEquation(0, cTReEqn);
        cTEqn.insertEquation(3, cTImEqn);

        // Coupling coeffs
        scalarField omegasf(cMesh.nCells(), omega.value());

        // Add coupling terms for TRe-TIm coupling
        cTEqn.insertEquationCoupling(0, 3,  omegasf);
        cTEqn.insertEquationCoupling(3, 0, -omegasf);
        cTEqn.insertEquationCoupling(1, 4,  omegasf);
        cTEqn.insertEquationCoupling(4, 1, -omegasf);
        cTEqn.insertEquationCoupling(2, 5,  omegasf);
        cTEqn.insertEquationCoupling(5, 2, -omegasf);

        // Block coupled solver call to solve T
        cTEqn.solve();

        // Retrieve solution for T
        cTEqn.retrieveSolution(0, cTRe.internalField());
        cTEqn.retrieveSolution(3, cTIm.internalField());

        cTRe.correctBoundaryConditions();
        cTIm.correctBoundaryConditions();

        // Calculate induced current
        cjRe == fvc::curl(cTRe);
        cjIm == fvc::curl(cTIm);

Info << "DEBUG (c): Reverse mapping" << endl;

        // Reverse map fields
        forAll (cMesh.C(), cellI)
        {
            SRe[cCellMap[cellI]] = cSRe[cellI];
            SIm[cCellMap[cellI]] = cSIm[cellI];

            TRe[cCellMap[cellI]] = cTRe[cellI];
            TIm[cCellMap[cellI]] = cTIm[cellI];

            jRe[cCellMap[cellI]] = cjRe[cellI];
            jIm[cCellMap[cellI]] = cjIm[cellI];
        }

        // debug
        cMesh.write();

        cSigma.write();
        cSRe.write();
        cSIm.write();
        cTRe.write();
        cTIm.write();
        cTsRe.write();
        cTsIm.write();
        cjRe.write();
        cjIm.write();
    }

//     {
//         fvScalarMatrix SReEqn
//         (
//             fvm::laplacian(SRe)
//           - fvc::div(TsRe)
//         );
//
//         fvScalarMatrix SImEqn
//         (
//             fvm::laplacian(SIm)
//           - fvc::div(TsIm)
//         );
//
//         // Set values for non-c cells
//         scalarField cCellSReValues(nCells.size(), 0.0); // TODO subfield
//         scalarField cCellSImValues(nCells.size(), 0.0); // TODO subfield
//         SReEqn.setValues(nCells, cCellSReValues);
//         SImEqn.setValues(nCells, cCellSImValues);
//
// //         // Set references
// //         label refCell = mesh.findCell(point(0.00, 0.00, 0.01));
// //         scalar refValue = 0.0;
// //
// //         SReEqn.setReference(refCell, refValue);
// //         SImEqn.setReference(refCell, refValue);
//
//         SReEqn.solve();
//         SImEqn.solve();
//
//         SRe.correctBoundaryConditions();
//         SIm.correctBoundaryConditions();
//     }

Info << "DEBUG: Solve S" << endl;

    {
        fvScalarMatrix SReEqn
        (
            fvm::laplacian(SRe)
          - fvc::div(TRe)
          - fvc::div(TsRe)
        );

        fvScalarMatrix SImEqn
        (
            fvm::laplacian(SIm)
          - fvc::div(TIm)
          - fvc::div(TsIm)
        );

//         // Set values for conducting cells
//         scalarField cCellSReValues(cCells.size(), 0.0);
//         scalarField cCellSImValues(cCells.size(), 0.0);
//         forAll (cCells, cellI)
//         {
//             cCellSReValues[cellI] = SRe.internalField()[cCells[cellI]];
//             cCellSImValues[cellI] = SIm.internalField()[cCells[cellI]];
//         }
//         SReEqn.setValues(cCells, cCellSReValues);
//         SImEqn.setValues(cCells, cCellSImValues);

        // Set references
        label refCell = mesh.findCell(point(0.00, 0.00, 0.01));
        scalar refValue = 0.0;

        SReEqn.setReference(refCell, refValue);
        SImEqn.setReference(refCell, refValue);

        SReEqn.solve();
        SImEqn.solve();

        SRe.correctBoundaryConditions();
        SIm.correctBoundaryConditions();
    }

//
//
//
//         fvVectorMatrix TReEqn
//         (
//             mgamma * fvm::laplacian(TRe)
// //           - (fvc::grad(mgamma) ^ fvc::curl(TRe))
// //           - mgamma * fvc::grad(fvc::div(TRe))
// //           + omega * TIm
//           + omega * TsIm
//           - omega * fvc::grad(SIm)
//         );
//
//         fvVectorMatrix TImEqn
//         (
//             mgamma * fvm::laplacian(TIm)
// //           - (fvc::grad(mgamma) ^ fvc::curl(TIm))
// //           - mgamma * fvc::grad(fvc::div(TIm))
// //           - omega * TRe
//           - omega * TsRe
//           + omega * fvc::grad(SRe)
//         );
//
//         // Set values for non-c cells
//         vectorField nCellTValues(nCells.size(), vector::zero);
//         TReEqn.setValues(nCells, nCellTValues);
//         TImEqn.setValues(nCells, nCellTValues);
//
//         // Prepare block system
//         fvBlockMatrix<vector6> TEqn(T);
//
//         // Insert equations for T
//         TEqn.insertEquation(0, TReEqn);
//         TEqn.insertEquation(3, TImEqn);
//
//         // Coupling coeffs
//         scalarField omegasf(mesh.nCells(), omega.value());
//
//         // Add coupling terms for TRe-TIm coupling
//         TEqn.insertEquationCoupling(0, 3,  omegasf);
//         TEqn.insertEquationCoupling(3, 0, -omegasf);
//         TEqn.insertEquationCoupling(1, 4,  omegasf);
//         TEqn.insertEquationCoupling(4, 1, -omegasf);
//         TEqn.insertEquationCoupling(2, 5,  omegasf);
//         TEqn.insertEquationCoupling(5, 2, -omegasf);
//
//         // Block coupled solver call to solve T
//         TEqn.solve();
//
//         // Retrieve solution
//         TEqn.retrieveSolution(0, TRe.internalField());
//         TEqn.retrieveSolution(3, TIm.internalField());
//
//         TRe.correctBoundaryConditions();
//         TIm.correctBoundaryConditions();
//
//
//
// //         fvScalarMatrix GReEqn
// //         (
// //             fvm::laplacian(GRe)
// //           - fvc::div(TRe)
// //         );
// //
// //         fvScalarMatrix GImEqn
// //         (
// //             fvm::laplacian(GIm)
// //           - fvc::div(TIm)
// //         );
// //
// //         // Set values for non-c cells
// //         scalarField nCellGValues(nCells.size(), 0.0);
// //         GReEqn.setValues(nCells, nCellGValues);
// //         GImEqn.setValues(nCells, nCellGValues);
// //
// //         GReEqn.solve();
// //         GImEqn.solve();
// //
// //         GRe.correctBoundaryConditions();
// //         GIm.correctBoundaryConditions();
// //
// //         TRe -= fvc::grad(GRe);
// //         TIm -= fvc::grad(GIm);
// //
// //         TRe.correctBoundaryConditions();
// //         TIm.correctBoundaryConditions();
//     }

    BRe == mu0 * (TRe + TsRe - fvc::grad(SRe));
    BIm == mu0 * (TIm + TsIm - fvc::grad(SIm));


