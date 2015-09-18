
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



    // Linear interpolate the magnetic scalar potential
    // to the faces (sS). We know that sS is obviously
    // wrong at the c/n-interface as the gradient of S
    // depends on the current vector potential (Ts) there.
    surfaceScalarField sSRe("sSRe", linearInterpolate(SRe));
    surfaceScalarField sSIm("sSIm", linearInterpolate(SIm));

    // This scope contains the correction procedure
    // of the face values of the interpolated magnetic
    // scalar sS at the c/n-interface. The whole code is
    // only valid for (wrongly applied) LINEAR interpolation.
    {
        // Mesh cell-face addressing
        const unallocLabelList& own = mesh.owner();
        const unallocLabelList& nei = mesh.neighbour();

        // Mesh geometry
        // (Sf/nf always points from owner to neighbour)
        const vectorField& Cf = mesh.Cf();
        const vectorField& C = mesh.C();
        const vectorField nf = mesh.Sf()/mesh.magSf();

        // References to internal fields of sS
        scalarField& sSReIn = sSRe.internalField();
        scalarField& sSImIn = sSIm.internalField();

        // Constant references to internal fields of T
        const vectorField& TReIn = TRe.internalField();
        const vectorField& TImIn = TIm.internalField();

        // Owner: conducting
        // Neighbour: non-conducting
        forAll (cnPNfaces.toc(), facei)
        {
            label faceI = cnPNfaces.toc()[facei];

            // Calculate Pf and fN
            scalar Pf = mag(nf[faceI] & (Cf[faceI] - C[own[faceI]]));
            scalar fN = mag(nf[faceI] & (C[nei[faceI]] - Cf[faceI]));

            // Calculate weight for T
            scalar wT = Pf*fN/(Pf+fN);

            // Correct face values of sS (plus!)
            // Note: Owner is in conducting region
            sSReIn[faceI] += wT * (nf[faceI] & TReIn[own[faceI]]);
            sSImIn[faceI] += wT * (nf[faceI] & TImIn[own[faceI]]);
        }

        // Owner: non-conducting
        // Neighbour: conducting
        forAll (cnNPfaces.toc(), facei)
        {
            label faceI = cnNPfaces.toc()[facei];

            // Calculate Pf and fN
            scalar Pf = mag(nf[faceI] & (Cf[faceI] - C[own[faceI]]));
            scalar fN = mag(nf[faceI] & (C[nei[faceI]] - Cf[faceI]));

            // Calculate weight for T
            scalar wT = Pf*fN/(Pf+fN);

            // Correct face values of sS (minus!)
            // Note: Neighbour is in conducting region
            sSReIn[faceI] -= wT * (nf[faceI] & TReIn[nei[faceI]]);
            sSImIn[faceI] -= wT * (nf[faceI] & TImIn[nei[faceI]]);
        }
    }

    // Calculate and store the gradient of the
    // corrected magnetic scalar potential (sS) on
    // the faces
    volVectorField SReGrad("SReGrad", fvc::grad(sSRe));
    volVectorField SImGrad("SImGrad", fvc::grad(sSIm));

    // TODO: To be removed after debug
    SReGrad.write();
    SImGrad.write();



    // This scope is to solve the current vector
    // potential (T) in the conducting region. Therefore
    // we will create a submesh and map all necessary
    // fields to this submesh. The boundary conditions
    // for T are very unusual but they can be realized
    // by means of the "newDirectionMixedFvPatchFields"
    // from foam-extend. To be more specific we need to set
    // the tangential components of T to "zero" (Txn = 0) and
    // apply a zero-gradient in normal direction (dT/dn = 0).
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



        // Create fields
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
                // in order to get zero tangential values
                // at the boundaries (T^n=0)
                cTReP.tangentialValueFraction() = 1.0;

                // Set normal value fraction to 0
                // to get zero normal derivative
                // which gives div(T)=0 with T^n=0
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
                // in order to get zero tangential values
                // at the boundaries (T^n=0)
                cTImP.tangentialValueFraction() = 1.0;

                // Set normal value fraction to 0
                // to get zero normal derivative
                // which gives div(T)=0 with T^n=0
                cTImP.normalValueFraction() = 0.0;
            }
        }



        // Map internal values of volume fields
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

        // Correct boundary values of volume fields
        cSigma.correctBoundaryConditions();

        cSReGrad.correctBoundaryConditions();
        cSImGrad.correctBoundaryConditions();

        cTRe.correctBoundaryConditions();
        cTIm.correctBoundaryConditions();

        cTsRe.correctBoundaryConditions();
        cTsIm.correctBoundaryConditions();


        // Real part
        fvVectorMatrix cTReEqn
        (
            rMu0/cSigma * fvm::laplacian(cTRe)
          - (fvc::grad(rMu0/cSigma) ^ fvc::curl(cTRe)) // TODO: Testing!
//           + omega * cTIm                            // Block-coupling
          + omega * cTsIm
          - omega * cSImGrad
        );

        // Imaginary part
        fvVectorMatrix cTImEqn
        (
            rMu0/cSigma * fvm::laplacian(cTIm)
          - (fvc::grad(rMu0/cSigma) ^ fvc::curl(cTIm)) // TODO: Testing!
//           - omega * cTRe                            // Block-coupling
          - omega * cTsRe
          + omega * cSReGrad
        );

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



        // Reverse map fields
        forAll (cMesh.C(), cellI)
        {
            TRe[cCellMap[cellI]] = cTRe[cellI];
            TIm[cCellMap[cellI]] = cTIm[cellI];

            jRe[cCellMap[cellI]] = cjRe[cellI];
            jIm[cCellMap[cellI]] = cjIm[cellI];
        }

        // TODO: Remove after debug!
        cMesh.write();
        cSigma.write();
        cSReGrad.write();
        cSImGrad.write();
        cTRe.write();
        cTIm.write();
        cTsRe.write();
        cTsIm.write();
        cjRe.write();
        cjIm.write();
    }



    // Calculate/update magnetic field density
    BRe == mu0 * (TRe + TsRe - SReGrad);
    BIm == mu0 * (TIm + TsIm - SImGrad);


