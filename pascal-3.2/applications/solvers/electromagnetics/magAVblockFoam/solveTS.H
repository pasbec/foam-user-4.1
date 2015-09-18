
Info << "DEBUG: Solve Ts" << endl;

    // Calculate source current vector potential (Ts)
    // from source given external source current
    // density (js). Ts is only depending on the mesh and
    // the current density distribution in the inductor
    // coils. As boundary conditions for Ts at "infinity"
    // we will use "slip" (normal magnetic)
    //
    // Note: Alternatively, the same could be realized by
    //       using the Biot-Savart-Law
    {
        // Real part
        fvVectorMatrix TsReEqn
        (
            fvm::laplacian(TsRe)
          + fvc::curl(jsRe)
        );

        // Imaginary part
        fvVectorMatrix TsImEqn
        (
            fvm::laplacian(TsIm)
          + fvc::curl(jsIm)
        );

        // Solve Ts equations
        TsReEqn.solve();
        TsImEqn.solve();

        TsRe.correctBoundaryConditions();
        TsIm.correctBoundaryConditions();
    }

Info << "DEBUG: Solve S" << endl;

    // Calculate magnetic scalar potential (S)
    // from source current vector potential and
    // current vector potential in the conducting region
    // from the last time step. As boundary conditions
    // at "infinity" we will use "zeroGradient" which is
    // consistent with "slip" for Ts.
    {
        // Real part
        fvScalarMatrix SReEqn
        (
            fvm::laplacian(SRe)
          - fvc::div(TRe)
          - fvc::div(TsRe)
        );

        // Imaginary part
        fvScalarMatrix SImEqn
        (
            fvm::laplacian(SIm)
          - fvc::div(TIm)
          - fvc::div(TsIm)
        );

        // Set references
        label refCell = mesh.findCell(point(0.00, 0.00, 0.01));
        scalar refValue = 0.0;

        SReEqn.setReference(refCell, refValue);
        SImEqn.setReference(refCell, refValue);

        // Solve S equations
        SReEqn.solve();
        SImEqn.solve();

        SRe.correctBoundaryConditions();
        SIm.correctBoundaryConditions();
    }

    // Calculate and store the gradient of
    // the magnetic scalar potential (S)
    volVectorField SReGrad("SReGrad", fvc::grad(SRe));
    volVectorField SImGrad("SImGrad", fvc::grad(SIm));

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

Info << "DEBUG: Create subMesh for conducting region" << endl;

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

Info << "DEBUG (c): Field creation for subMesh" << endl;

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
//             zeroGradientFvPatchField<scalar>::typeName
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
//             zeroGradientFvPatchField<scalar>::typeName
//         );

        volVectorField cSReGrad
        (
            IOobject
            (
                "SReGrad",
                cMesh.time().timeName(),
                cMesh,
                IOobject::NO_READ
            ),
            cMesh,
            dimensioned<vector>(word(), SReGrad.dimensions(), pTraits<vector>::zero),
            zeroGradientFvPatchField<vector>::typeName
        );

        volVectorField cSImGrad
        (
            IOobject
            (
                "SImGrad",
                cMesh.time().timeName(),
                cMesh,
                IOobject::NO_READ
            ),
            cMesh,
            dimensioned<vector>(word(), SImGrad.dimensions(), pTraits<vector>::zero),
            zeroGradientFvPatchField<vector>::typeName
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

Info << "DEBUG (c): Boundary conditions for subMesh fields" << endl;

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

Info << "DEBUG (c): Mapping to subMesh" << endl;

        // Map volume fields internal
        forAll (cMesh.C(), cellI)
        {
            cSigma.internalField()[cellI] = sigma.internalField()[cCellMap[cellI]];

            cSReGrad.internalField()[cellI] = SReGrad.internalField()[cCellMap[cellI]];
            cSImGrad.internalField()[cellI] = SImGrad.internalField()[cCellMap[cellI]];

            cTRe.internalField()[cellI] = TRe.internalField()[cCellMap[cellI]];
            cTIm.internalField()[cellI] = TIm.internalField()[cCellMap[cellI]];

            cTsRe.internalField()[cellI] = TsRe.internalField()[cCellMap[cellI]];
            cTsIm.internalField()[cellI] = TsIm.internalField()[cCellMap[cellI]];
        }

        // Correct boundary conditions
        cSigma.correctBoundaryConditions();

        cSReGrad.correctBoundaryConditions();
        cSImGrad.correctBoundaryConditions();

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

Info << "DEBUG (c): Write subMesh fields with explicit influence to T equation" << endl;

        cMesh.write();

        cSigma.write();

        cTsRe.write();
        cTsIm.write();

        cSReGrad.write();
        cSImGrad.write();

Info << "DEBUG (c): Prepare matrices for T" << endl;

        fvVectorMatrix cTReEqn
        (
            rMu0/cSigma * fvm::laplacian(cTRe)
          - (fvc::grad(rMu0/cSigma) ^ fvc::curl(cTRe)) // TODO: Testing!
//           + omega * cTIm
          + omega * cTsIm
          - omega * cSImGrad
        );

        fvVectorMatrix cTImEqn
        (
            rMu0/cSigma * fvm::laplacian(cTIm)
          - (fvc::grad(rMu0/cSigma) ^ fvc::curl(cTIm)) // TODO: Testing!
//           - omega * cTRe
          - omega * cTsRe
          + omega * cSReGrad
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

Info << "DEBUG (c): Reverse mapping from subMesh" << endl;

        // Reverse map fields
        forAll (cMesh.C(), cellI)
        {
            TRe[cCellMap[cellI]] = cTRe[cellI];
            TIm[cCellMap[cellI]] = cTIm[cellI];

            jRe[cCellMap[cellI]] = cjRe[cellI];
            jIm[cCellMap[cellI]] = cjIm[cellI];
        }

Info << "DEBUG (c): Write subMesh fields T and dependencies" << endl;

        cTRe.write();
        cTIm.write();
        cjRe.write();
        cjIm.write();
    }

    BRe == mu0 * (TRe + TsRe - SReGrad);
    BIm == mu0 * (TIm + TsIm - SImGrad);


