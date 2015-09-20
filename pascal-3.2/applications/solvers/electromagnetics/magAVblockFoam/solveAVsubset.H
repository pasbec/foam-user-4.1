
    {
        // Prepare A-equations
        fvVectorMatrix AReEqnSubset = AReEqn;
        fvVectorMatrix AImEqnSubset = AImEqn;

        // Add explicit parts of G
        AReEqnSubset += alpha * fvc::grad(GIm);
        AImEqnSubset -= alpha * fvc::grad(GRe);

        // Add explicit parts of V
        AReEqnSubset -= sigma * fvc::grad(VRe);
        AImEqnSubset -= sigma * fvc::grad(VIm);

        // Prepare block system
        fvBlockMatrix<vector6> AEqn(A);

        // Insert equations for A
        AEqn.insertEquation(0, AReEqnSubset);
        AEqn.insertEquation(3, AImEqnSubset);

        // References to internal alpha field
        const scalarField& alphaIn = alpha.internalField();

        // Add coupling terms for ARe-AIm coupling
        AEqn.insertEquationCoupling(0, 3,  alphaIn);
        AEqn.insertEquationCoupling(3, 0, -alphaIn);
        AEqn.insertEquationCoupling(1, 4,  alphaIn);
        AEqn.insertEquationCoupling(4, 1, -alphaIn);
        AEqn.insertEquationCoupling(2, 5,  alphaIn);
        AEqn.insertEquationCoupling(5, 2, -alphaIn);

        // Block coupled solver call to solve A
        AEqn.solve();

        // Retrieve solution
        AEqn.retrieveSolution(0, ARe.internalField());
        AEqn.retrieveSolution(3, AIm.internalField());

        ARe.correctBoundaryConditions();
        AIm.correctBoundaryConditions();

        // Create interpolated surface normal components of A
        surfaceScalarField snARe("snARe", fvc::interpolate(ARe) & mesh.Sf()/mesh.magSf());
        surfaceScalarField snAIm("snAIm", fvc::interpolate(AIm) & mesh.Sf()/mesh.magSf());

        //
        fvScalarMatrix GReEqn(fvm::laplacian(GRe) + fvc::div(ARe));
        fvScalarMatrix GImEqn(fvm::laplacian(GIm) + fvc::div(AIm));

        // Solve V
        GReEqn.solve();
        GImEqn.solve();

        GRe.correctBoundaryConditions();
        GIm.correctBoundaryConditions();

        // Create interpolated surface gradient of G
        surfaceScalarField sGradGRe
        (
            "sGradGRe",
            fvc::interpolate(fvc::grad(GRe)) & mesh.Sf()
        );

        surfaceScalarField sGradGIm
        (
            "sGradGIm",
            fvc::interpolate(fvc::grad(GIm)) & mesh.Sf()
        );

        // Create set for conductor cells
        labelHashSet conductorCellsHashSet;
        forAll(mesh.cells(), cellI)
        {
            if (sigma.internalField()[cellI] > SMALL)
            {
                conductorCellsHashSet.insert(cellI);
            }
        }

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
        cMeshSubSet.setLargeCellSubset(conductorCellsHashSet);

        // Create sub-mesh
        fvMesh& cMesh = cMeshSubSet.subMesh();

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
            dimensionedVector(word(), ARe.dimensions(), vector::zero)
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
            dimensionedVector(word(), AIm.dimensions(), vector::zero)
        );

        volScalarField cGRe
        (
            IOobject
            (
                "GRe",
                cMesh.time().timeName(),
                cMesh,
                IOobject::NO_READ
            ),
            cMesh,
            dimensionedScalar(word(), GRe.dimensions(), 0.0),
            fixedValueFvPatchField<scalar>::typeName
        );

        volScalarField cGIm
        (
            IOobject
            (
                "GIm",
                cMesh.time().timeName(),
                cMesh,
                IOobject::NO_READ
            ),
            cMesh,
            dimensionedScalar(word(), GIm.dimensions(), 0.0),
            fixedValueFvPatchField<scalar>::typeName
        );

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
            dimensionedScalar(word(), VRe.dimensions(), 0.0),
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
            dimensionedScalar(word(), VIm.dimensions(), 0.0),
            fixedGradientFvPatchField<scalar>::typeName
        );

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
            dimensionedScalar(word(), sigma.dimensions(), 0.0)
        );

        volVectorField cGradAlpha
        (
            IOobject
            (
                "gradAlpha",
                cMesh.time().timeName(),
                cMesh,
                IOobject::NO_READ
            ),
            cMesh,
            dimensionedVector(word(), gradAlpha.dimensions(), vector::zero)
        );

        // Map volume fields
        forAll (cMesh.C(), cellI)
        {
            cARe.internalField()[cellI] = ARe.internalField()[cMeshSubSet.cellMap()[cellI]];
            cAIm.internalField()[cellI] = AIm.internalField()[cMeshSubSet.cellMap()[cellI]];

            cGRe.internalField()[cellI] = GRe.internalField()[cMeshSubSet.cellMap()[cellI]];
            cGIm.internalField()[cellI] = GIm.internalField()[cMeshSubSet.cellMap()[cellI]];

            cVRe.internalField()[cellI] = VRe.internalField()[cMeshSubSet.cellMap()[cellI]];
            cVIm.internalField()[cellI] = VIm.internalField()[cMeshSubSet.cellMap()[cellI]];

            cSigma.internalField()[cellI] = sigma.internalField()[cMeshSubSet.cellMap()[cellI]];
            cGradAlpha.internalField()[cellI] = gradAlpha.internalField()[cMeshSubSet.cellMap()[cellI]];
        }

        surfaceScalarField cSnARe
        (
            IOobject
            (
                "snARe",
                cMesh.time().timeName(),
                cMesh,
                IOobject::NO_READ
            ),
            cMesh,
            dimensionedScalar(word(), snARe.dimensions(), 0.0)
        );

        surfaceScalarField cSnAIm
        (
            IOobject
            (
                "snAIm",
                cMesh.time().timeName(),
                cMesh,
                IOobject::NO_READ
            ),
            cMesh,
            dimensionedScalar(word(), snAIm.dimensions(), 0.0)
        );

        surfaceScalarField cSgradGRe
        (
            IOobject
            (
                "sGradGRe",
                cMesh.time().timeName(),
                cMesh,
                IOobject::NO_READ
            ),
            cMesh,
            dimensionedScalar(word(), sGradGRe.dimensions(), 0.0)
        );

        surfaceScalarField cSgradGIm
        (
            IOobject
            (
                "sGradGIm",
                cMesh.time().timeName(),
                cMesh,
                IOobject::NO_READ
            ),
            cMesh,
            dimensionedScalar(word(), sGradGIm.dimensions(), 0.0)
        );

        // Map surface fields
        forAll (cMesh.owner(), faceI)
        {
            cSnARe.internalField()[faceI] = snARe.internalField()[cMeshSubSet.faceMap()[faceI]];
            cSnAIm.internalField()[faceI] = snAIm.internalField()[cMeshSubSet.faceMap()[faceI]];

            cSgradGRe.internalField()[faceI] = sGradGRe.internalField()[cMeshSubSet.faceMap()[faceI]];
            cSgradGIm.internalField()[faceI] = sGradGIm.internalField()[cMeshSubSet.faceMap()[faceI]];
        }
        forAll(cMesh.boundaryMesh(), patchI)
        {
            polyPatch patch = cMesh.boundaryMesh()[patchI];

            forAll(patch, facei)
            {
                label faceI = patch.start() + facei;

                // FIXME: This does only work, if boundary faces from
                //        subsetmesh do NOT coincide with b-faces of base mesh
                cSnARe.boundaryField()[patchI][facei] = snARe[cMeshSubSet.faceMap()[faceI]];
                cSnAIm.boundaryField()[patchI][facei] = snAIm[cMeshSubSet.faceMap()[faceI]];

                cSgradGRe.boundaryField()[patchI][facei] = sGradGRe[cMeshSubSet.faceMap()[faceI]];
                cSgradGIm.boundaryField()[patchI][facei] = sGradGIm[cMeshSubSet.faceMap()[faceI]];
            }
        }

        // Boundary conditions for V
        forAll(cVRe.boundaryField(), patchI)
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
                cVReP.gradient() =   omega.value() * cSnAIm.boundaryField()[patchI]
                                   + omega.value() * cSgradGIm.boundaryField()[patchI] / cMesh.magSf().boundaryField()[patchI];
            }
        }
        forAll(cVIm.boundaryField(), patchI)
        {
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
                cVImP.gradient() = - omega.value() * cSnARe.boundaryField()[patchI]
                                   - omega.value() * cSgradGRe.boundaryField()[patchI] / cMesh.magSf().boundaryField()[patchI];
            }
        }

        // Interpolate sigma to face centers
        surfaceScalarField cSigmaf
        (
            "sigmaf",
            fvc::interpolate(cSigma, "interpolate(sigma)")
        );

        // Prepare V-equations
        fvScalarMatrix cVReEqn(fvm::laplacian(cSigmaf, cVRe, "laplacian(VRe)"));
        fvScalarMatrix cVImEqn(fvm::laplacian(cSigmaf, cVIm, "laplacian(VIm)"));

        surfaceScalarField sigmaGradGRe("sigmaGradGRe", cSigmaf*cSgradGRe);
        surfaceScalarField sigmaGradGIm("sigmaGradGIm", cSigmaf*cSgradGIm);

        cVReEqn -= omega*fvc::div(sigmaGradGIm);
        cVImEqn += omega*fvc::div(sigmaGradGRe);

//         cVReEqn -= omega*fvc::laplacian(cSigmaf, cGIm, "laplacian(GIm)");
//         cVImEqn += omega*fvc::laplacian(cSigmaf, cGRe, "laplacian(GRe)");

        cVReEqn -= cGradAlpha&cAIm;
        cVImEqn += cGradAlpha&cARe;

        // Set references
        label refCell = cMesh.findCell(point(0.00, 0.00, 0.01));
        scalar refValue = 0.0;

        cVReEqn.setReference(refCell, refValue);
        cVImEqn.setReference(refCell, refValue);

        // Solve V-equations
        cVReEqn.solve();
        cVImEqn.solve();

        cVRe -= average(cVRe);
        cVIm -= average(cVIm);

        // Reverse map V
        forAll (cVRe.internalField(), cellI)
        {
            VRe[cMeshSubSet.cellMap()[cellI]] = cVRe[cellI];
        }
        forAll (cVIm.internalField(), cellI)
        {
            VIm[cMeshSubSet.cellMap()[cellI]] = cVIm[cellI];
        }

        cMesh.write();
        cVRe.write();
        cVIm.write();
        cSnARe.write();
        cSnAIm.write();
    }

    jRe ==   alpha * AIm + alpha * fvc::grad(GIm) - sigma * fvc::grad(VRe);
    jIm == - alpha * ARe - alpha * fvc::grad(GRe) - sigma * fvc::grad(VIm);

    BRe == fvc::curl(ARe);
    BIm == fvc::curl(AIm);
