
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

    // Create set for non-conducting cells
    labelHashSet nonConductingCellsHashSet;
    forAll(mesh.cells(), cellI)
    {
        if (sigma.internalField()[cellI] < SMALL)
        {
            nonConductingCellsHashSet.insert(cellI);
        }
    }
    labelList nonConductingCells = nonConductingCellsHashSet.toc();

    // Create set for conducting cells
    labelHashSet conductingCellsHashSet;
    forAll(mesh.cells(), cellI)
    {
        if (sigma.internalField()[cellI] > SMALL)
        {
            conductingCellsHashSet.insert(cellI);
        }
    }
    labelList conductingCells = conductingCellsHashSet.toc();

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
        cMeshSubSet.setLargeCellSubset(conductingCellsHashSet);

        // Create sub-mesh
        fvMesh& cMesh = cMeshSubSet.subMesh();

Info << "DEBUG: Fields" << endl;

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
            calculatedFvPatchField<scalar>::typeName
        );

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
            fixedGradientFvPatchField<scalar>::typeName
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
            fixedGradientFvPatchField<scalar>::typeName
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
            calculatedFvPatchField<vector>::typeName
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
            calculatedFvPatchField<vector>::typeName
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

Info << "DEBUG: Boundary conditions" << endl;

        // Boundary conditions
        forAll(cMesh.boundary(), patchI)
        {
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

        cTRe.write();
        cTIm.write();

Info << "DEBUG: Mapping" << endl;

        // Map volume fields
        forAll (cMesh.C(), cellI)
        {
            cSigma.internalField()[cellI] = sigma.internalField()[cMeshSubSet.cellMap()[cellI]];

            cSRe.internalField()[cellI] = SRe.internalField()[cMeshSubSet.cellMap()[cellI]];
            cSIm.internalField()[cellI] = SIm.internalField()[cMeshSubSet.cellMap()[cellI]];

            cTRe.internalField()[cellI] = TRe.internalField()[cMeshSubSet.cellMap()[cellI]];
            cTIm.internalField()[cellI] = TIm.internalField()[cMeshSubSet.cellMap()[cellI]];

            cTsRe.internalField()[cellI] = TsRe.internalField()[cMeshSubSet.cellMap()[cellI]];
            cTsIm.internalField()[cellI] = TsIm.internalField()[cMeshSubSet.cellMap()[cellI]];
        }
        forAll(cMesh.boundaryMesh(), patchI)
        {
            polyPatch patch = cMesh.boundaryMesh()[patchI];

            forAll(patch, facei)
            {
                label faceI = patch.start() + facei;
            }

            // Use internalField for boundary field for sigma
            cSigma.boundaryField()[patchI] = cSigma.boundaryField()[patchI].patchInternalField();
        }

Info << "DEBUG: Marices" << endl;

        fvVectorMatrix cTReEqn
        (
            rMu0/cSigma * fvm::laplacian(cTRe)
//           + omega * cTIm
          + omega * cTsIm
//           - omega * fvc::grad(cSIm)
        );

        fvVectorMatrix cTImEqn
        (
            rMu0/cSigma * fvm::laplacian(cTIm)
//           - omega * cTRe
          - omega * cTsRe
//           + omega * fvc::grad(cSRe)
        );

Info << "DEBUG: Blockmatrix" << endl;

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

        // Retrieve solution
        cTEqn.retrieveSolution(0, cTRe.internalField());
        cTEqn.retrieveSolution(3, cTIm.internalField());

        cTRe.correctBoundaryConditions();
        cTIm.correctBoundaryConditions();

        cjRe == fvc::curl(cTRe);
        cjIm == fvc::curl(cTIm);

Info << "DEBUG Reverse mapping" << endl;

        // Reverse map T
        forAll (cMesh.C(), cellI)
        {
            TRe[cMeshSubSet.cellMap()[cellI]] = cTRe[cellI];
            TIm[cMeshSubSet.cellMap()[cellI]] = cTIm[cellI];

            jRe[cMeshSubSet.cellMap()[cellI]] = cjRe[cellI];
            jIm[cMeshSubSet.cellMap()[cellI]] = cjIm[cellI];
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
//           - fvc::div(TRe)
//           - fvc::div(TsRe)
//         );
//
//         fvScalarMatrix SImEqn
//         (
//             fvm::laplacian(SIm)
//           - fvc::div(TIm)
//           - fvc::div(TsIm)
//         );
//
// //         // Set values for non-conducting cells
// //         scalarField nonConductingCellSValues(nonConductingCells.size(), 0.0);
// //         SReEqn.setValues(nonConductingCells, nonConductingCellSValues);
// //         SImEqn.setValues(nonConductingCells, nonConductingCellSValues);
//
//         // Set references
//         label refCell = mesh.findCell(point(0.00, 0.00, 0.01));
//         scalar refValue = 0.0;
//
//         SReEqn.setReference(refCell, refValue);
//         SImEqn.setReference(refCell, refValue);
//
//         SReEqn.solve();
//         SImEqn.solve();
//
//         SRe.correctBoundaryConditions();
//         SIm.correctBoundaryConditions();
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
//         // Set values for non-conducting cells
//         vectorField nonConductingCellTValues(nonConductingCells.size(), vector::zero);
//         TReEqn.setValues(nonConductingCells, nonConductingCellTValues);
//         TImEqn.setValues(nonConductingCells, nonConductingCellTValues);
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
// //         // Set values for non-conducting cells
// //         scalarField nonConductingCellGValues(nonConductingCells.size(), 0.0);
// //         GReEqn.setValues(nonConductingCells, nonConductingCellGValues);
// //         GImEqn.setValues(nonConductingCells, nonConductingCellGValues);
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


