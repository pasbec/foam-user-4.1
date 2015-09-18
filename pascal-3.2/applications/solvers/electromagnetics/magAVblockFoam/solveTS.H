
Info << "DEBUG: Define c/n cellSets and faceSets" << endl;

    // Create set for conducting cells
    labelHashSet cCellsHashSet;

    forAll(mesh.cells(), cellI)
    {
        if (sigma.internalField()[cellI] > SMALL)
        {
            cCellsHashSet.insert(cellI);
        }
    }

    cellSet cCells
    (
        mesh,
        "cellSet_solver_c",
        cCellsHashSet
    ); cCells.write();

    // Create set for non-conducting cells
    labelHashSet nCellsHashSet;

    forAll(mesh.cells(), cellI)
    {
        if (sigma.internalField()[cellI] < SMALL)
        {
            nCellsHashSet.insert(cellI);
        }
    }

    cellSet nCells
    (
        mesh,
        "cellSet_solver_n",
        nCellsHashSet
    ); nCells.write();

    // Create sets for conducting/non-conducting faces (cn):
    // cnPN) Owner is in conducting region and neighbour
    //       lies in non-conducting region
    // cnNP) Owner is in non-conducting region and
    //       neighbour lies in conducting region
    labelHashSet cnFacesHashSet;
    labelHashSet cnPNfacesHashSet;
    labelHashSet cnNPfacesHashSet;

    {
        // Mesh references
        const unallocLabelList& own = mesh.owner();
        const unallocLabelList& nei = mesh.neighbour();

        // Find and insert faces
        forAll (own, faceI)
        {
            if
            (
                (sigma.internalField()[own[faceI]] > SMALL)
             && (sigma.internalField()[nei[faceI]] < SMALL)
            )
            {
                cnFacesHashSet.insert(faceI);
                cnPNfacesHashSet.insert(faceI);
            }

            if
            (
                (sigma.internalField()[own[faceI]] < SMALL)
             && (sigma.internalField()[nei[faceI]] > SMALL)
            )
            {
                cnFacesHashSet.insert(faceI);
                cnNPfacesHashSet.insert(faceI);
            }
        }
    }

    faceSet cnFaces
    (
        mesh,
        "cellSet_solver_cn",
        cnFacesHashSet
    ); cnFaces.write();

    faceSet cnPNfaces
    (
        mesh,
        "cellSet_solver_cn_PN",
        cnPNfacesHashSet
    ); cnPNfaces.write();

    faceSet cnNPfaces
    (
        mesh,
        "cellSet_solver_cn_NP",
        cnNPfacesHashSet
    ); cnNPfaces.write();

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
    // from source current vector potential (Ts) and
    // current vector potential (T) in the conducting region
    // from the last time step. As boundary conditions
    // at "infinity" we will use "zeroGradient" which is
    // consistent with "slip" for Ts as H = Ts - grad(S)
    // in order to get H=0 there.
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

Info << "DEBUG: Corrected gradient of S" << endl;

    // Linear interpolate the magnetic scalar potential
    // tho the faces (sS). We know that sS is completely
    // wrong at the c/n-interface as the gradient of S
    // depends on the current vector potential (Ts) there.
    surfaceScalarField sSRe("sSRe", linearInterpolate(SRe));
    surfaceScalarField sSIm("sSIm", linearInterpolate(SIm));

    // Correction of the face values of sS at the
    // c/n-interface
    {
        // Mesh references
        const unallocLabelList& own = mesh.owner();
        const unallocLabelList& nei = mesh.neighbour();

        // Mesh geometry references
        // (Sf always points from owner to neighbour)
        const vectorField& Cf = mesh.Cf();
        const vectorField& C = mesh.C();
        const vectorField& Sf = mesh.Sf();
        const scalarField& magSf = mesh.magSf();

        // References to internal fields of sS
        scalarField& sSReIn = sSRe.internalField();
        scalarField& sSImIn = sSIm.internalField();

        // References to internal fields of sS
        const vectorField& TReIn = TRe.internalField();
        const vectorField& TImIn = TIm.internalField();

        // Owner: conducting
        // Neighbour: non-conducting
        forAll (cnPNfaces.toc(), facei)
        {
            label faceI = cnPNfaces.toc()[facei];

            // Note: Copy from surfaceInterpolation.C
            // Calculate Pf (SfdOwn) and fN (SfdNei)
            scalar SfdOwn = mag(Sf[faceI] & (Cf[faceI] - C[own[faceI]]));
            scalar SfdNei = mag(Sf[faceI] & (C[nei[faceI]] - Cf[faceI]));

            // Calculate weight for T
            scalar wT = SfdOwn*SfdNei/(SfdOwn+SfdNei);

            // Grab normal part of T from owner (conducting region)
            scalar SfTRe = (Sf[faceI]/magSf[faceI]) & TReIn[own[faceI]];
            scalar SfTIm = (Sf[faceI]/magSf[faceI]) & TImIn[own[faceI]];

            // TODO Naming PN/NP should be cn/nc
// Info << "PN: facei = " << facei << endl;
// Info << "PN: faceI = " << faceI << endl;
// Info << "PN: own[faceI] = " << own[faceI] << endl;
// Info << "PN: nei[faceI] = " << nei[faceI] << endl;
// Info << "PN: wT = " << wT << endl;
// Info << "PN: gMax(TReIn) = " << gMax(TReIn) << endl;
// Info << "PN: gMax(TImIn) = " << gMax(TImIn) << endl;
// Info << "PN: gMin(TReIn) = " << gMin(TReIn) << endl;
// Info << "PN: gMin(TImIn) = " << gMin(TImIn) << endl;
// Info << "PN: TReIn[own[faceI]] = " << TReIn[own[faceI]] << endl;
// Info << "PN: TImIn[own[faceI]] = " << TImIn[own[faceI]] << endl;
// Info << "PN: TReIn[own[faceI]] = " << TReIn[nei[faceI]] << endl;
// Info << "PN: TImIn[own[faceI]] = " << TImIn[nei[faceI]] << endl;
// Info << "PN: SfTRe = " << SfTRe << endl;
// Info << "PN: SfTIm = " << SfTIm << endl;
// Info << "PN: sSReIn[faceI] = " << sSReIn[faceI] << endl;
// Info << "PN: sSImIn[faceI] = " << sSImIn[faceI] << endl;

            // Correct face values of sS (plus!)
            sSReIn[faceI] += wT * SfTRe;
            sSImIn[faceI] += wT * SfTIm;

// Info << "PN: sSReInCorr[faceI] = " << sSReIn[faceI] << endl;
// Info << "PN: sSImInCorr[faceI] = " << sSImIn[faceI] << endl;
        }

        // Owner: non-conducting
        // Neighbour: conducting
        forAll (cnNPfaces.toc(), facei)
        {
            label faceI = cnNPfaces.toc()[facei];

            // Copy from surfaceInterpolation.C
            // Calculate Pf (SfdOwn) and fN (SfdNei)
            scalar SfdOwn = mag(Sf[faceI] & (Cf[faceI] - C[own[faceI]]));
            scalar SfdNei = mag(Sf[faceI] & (C[nei[faceI]] - Cf[faceI]));

            // Calculate weight for T
            scalar wT = SfdOwn*SfdNei/(SfdOwn+SfdNei);

            // Grab normal part of T from neighbour (conducting region)
            scalar SfTRe = (Sf[faceI]/magSf[faceI]) & TReIn[nei[faceI]];
            scalar SfTIm = (Sf[faceI]/magSf[faceI]) & TImIn[nei[faceI]];

            // TODO Naming PN/NP should be cn/nc
// Info << "NP: facei = " << facei << endl;
// Info << "NP: faceI = " << faceI << endl;
// Info << "NP: own[faceI] = " << own[faceI] << endl;
// Info << "NP: nei[faceI] = " << nei[faceI] << endl;
// Info << "NP: wT = " << wT << endl;
// Info << "NP: gMax(TReIn) = " << gMax(TReIn) << endl;
// Info << "NP: gMax(TImIn) = " << gMax(TImIn) << endl;
// Info << "NP: gMin(TReIn) = " << gMin(TReIn) << endl;
// Info << "NP: gMin(TImIn) = " << gMin(TImIn) << endl;
// Info << "NP: TReIn[own[faceI]] = " << TReIn[own[faceI]] << endl;
// Info << "NP: TImIn[own[faceI]] = " << TImIn[own[faceI]] << endl;
// Info << "NP: TReIn[own[faceI]] = " << TReIn[nei[faceI]] << endl;
// Info << "NP: TImIn[own[faceI]] = " << TImIn[nei[faceI]] << endl;
// Info << "NP: SfTRe = " << SfTRe << endl;
// Info << "NP: SfTIm = " << SfTIm << endl;
// Info << "NP: sSReIn[faceI] = " << sSReIn[faceI] << endl;
// Info << "NP: sSImIn[faceI] = " << sSImIn[faceI] << endl;

            // Correct face values of sS (minus!)
            sSReIn[faceI] -= wT * SfTRe;
            sSImIn[faceI] -= wT * SfTIm;

// Info << "NP: sSReInCorr[faceI] = " << sSReIn[faceI] << endl;
// Info << "NP: sSImInCorr[faceI] = " << sSImIn[faceI] << endl;
        }
    }

    // Calculate and store the gradient of
    // the corrected magnetic scalar potential (sS)
    volVectorField SReGrad("SReGrad", fvc::grad(sSRe));
    volVectorField SImGrad("SImGrad", fvc::grad(sSIm));

Info << "DEBUG: Write gradient of sS" << endl;

    SReGrad.write();
    SImGrad.write();

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


