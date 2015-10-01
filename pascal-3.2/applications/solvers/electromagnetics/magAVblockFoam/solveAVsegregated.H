
    // Solve A
    {
        // Add explicit parts of electric scalar potential V
        ARexEqn -= sigma * VReGrad.component(0);
        AReyEqn -= sigma * VReGrad.component(1);
        ARezEqn -= sigma * VReGrad.component(2);
        AImxEqn -= sigma * VImGrad.component(0);
        AImyEqn -= sigma * VImGrad.component(1);
        AImzEqn -= sigma * VImGrad.component(2);

        // References to internal alpha field
        const scalarField& alphaIn = alpha.internalField();

        // Prepare block system for A
        fvBlockMatrix<vector6> AEqn(A);

        // Insert equations for A
        AEqn.insertEquation(0, ARexEqn);
        AEqn.insertEquation(1, AReyEqn);
        AEqn.insertEquation(2, ARezEqn);
        AEqn.insertEquation(3, AImxEqn);
        AEqn.insertEquation(4, AImyEqn);
        AEqn.insertEquation(5, AImzEqn);

        // Add coupling terms for ARe-AIm coupling
        AEqn.insertEquationCoupling(0, 3,  alphaIn);
        AEqn.insertEquationCoupling(3, 0, -alphaIn);
        AEqn.insertEquationCoupling(1, 4,  alphaIn);
        AEqn.insertEquationCoupling(4, 1, -alphaIn);
        AEqn.insertEquationCoupling(2, 5,  alphaIn);
        AEqn.insertEquationCoupling(5, 2, -alphaIn);

        // Solve A
        AEqn.solve();

        // Retrieve solution
        AEqn.retrieveSolution(0, ARex.internalField());
        AEqn.retrieveSolution(1, ARey.internalField());
        AEqn.retrieveSolution(2, ARez.internalField());
        AEqn.retrieveSolution(3, AImx.internalField());
        AEqn.retrieveSolution(4, AImy.internalField());
        AEqn.retrieveSolution(5, AImz.internalField());

        ARex.correctBoundaryConditions();
        ARey.correctBoundaryConditions();
        ARez.correctBoundaryConditions();
        AImx.correctBoundaryConditions();
        AImy.correctBoundaryConditions();
        AImz.correctBoundaryConditions();

        // Replace components
        ARe.replace(0,ARex);
        ARe.replace(1,ARey);
        ARe.replace(2,ARez);
        AIm.replace(0,AImx);
        AIm.replace(1,AImy);
        AIm.replace(2,AImz);

        ARe.correctBoundaryConditions();
        AIm.correctBoundaryConditions();

        // TODO: Remove after debug
        {
            vector6Field AResIn = AEqn.residual(A.internalField());

            volVectorField AReRes
            (
                IOobject
                (
                    "AReRes",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensioned<vector>(word(),dimless,pTraits<vector>::zero),
                fixedValueFvPatchField<vector>::typeName
            );
            volVectorField AImRes
            (
                IOobject
                (
                    "AImRes",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensioned<vector>(word(),dimless,pTraits<vector>::zero),
                fixedValueFvPatchField<vector>::typeName
            );

            AReRes.internalField().replace(0, AResIn.component(0));
            AReRes.internalField().replace(1, AResIn.component(1));
            AReRes.internalField().replace(2, AResIn.component(2));
            AImRes.internalField().replace(0, AResIn.component(4));
            AImRes.internalField().replace(1, AResIn.component(5));
            AImRes.internalField().replace(2, AResIn.component(6));

            AReRes.write();
            AImRes.write();
        }
    }

    // Solve V
    {
        // Interpolate A to face centers
        surfaceScalarField sARe
        (
            "sARe",
            fvc::interpolate(ARe) & mesh.Sf()
        );
        surfaceScalarField sAIm
        (
            "sAIm",
            fvc::interpolate(AIm) & mesh.Sf()
        );

        // Enforce zero normal components (ARe*n=0)
        // and (AIm*n=0) at the conductor interface as an
        // inner constraint
        {
            // Replace face values with zeros
            forAll (cnFaces, facei)
            {
                label faceI = cnFaces[facei];
                sARe[faceI] = 0.0;
                sAIm[faceI] = 0.0;
            }
        }

        // Coupling terms of magnetic vector potential A
        // in equation for electric scalar potential V
        {

//             VReEqn -= omega * fvc::div(sSigmaV*sAIm);
//             VImEqn += omega * fvc::div(sSigmaV*sARe);

            VReEqn -= omega * (fvc::grad(sSigmaV) & AIm);
            VImEqn += omega * (fvc::grad(sSigmaV) & ARe);

            VReEqn -= omega * sigmaV * fvc::div(sAIm);
            VImEqn += omega * sigmaV * fvc::div(sARe);
//             VReEqn -= omega * sigmaV * AImDiv;
//             VImEqn += omega * sigmaV * AReDiv;

            // Remove all source values from non-conducting cells
            VReEqn.source() *= sigmaIc;
            VImEqn.source() *= sigmaIc;
        }

        // Finish assembling of V
        {
            // Reference settings
            // TODO: Automatically search refcell in
            //       conducting region
            label refCell = mesh.findCell(point(0.00, 0.00, 0.01));
            scalar refValue = 0.0;

            // Force reference values for V-matrices in conducting
            // region.
            VReEqn.source()[refCell] +=  VReEqn.diag()[refCell]*refValue;
            VReEqn.diag()[refCell] +=  VReEqn.diag()[refCell];
            VImEqn.source()[refCell] +=  VImEqn.diag()[refCell]*refValue;
            VImEqn.diag()[refCell] +=  VImEqn.diag()[refCell];

            // Remove non-conducting cells from solution
            // TODO: In how far will this change upper/lower of
            // cells in contact with non-conducting region?
            // NOTE: The usage of setValues does not work
            //       correctly with block-matrices
//             scalarField nCellsV(nCells.size(), 0.0);
//             VReEqn.setValues(nCells, nCellsV);
//             VImEqn.setValues(nCells, nCellsV);
        }

        // Prepare block system for V
        fvBlockMatrix<vector2> VEqn(V);

        // Insert equations for V
        VEqn.insertEquation(0, VReEqn);
        VEqn.insertEquation(1, VImEqn);

        // Solve V
        VEqn.solve();

        // Retrieve solution
        VEqn.retrieveSolution(0, VRe.internalField());
        VEqn.retrieveSolution(1, VIm.internalField());

        VRe.correctBoundaryConditions();
        VIm.correctBoundaryConditions();

        // TODO: Remove after debug
        {
            vector2Field VResIn = VEqn.residual(V.internalField());

            volScalarField VReRes
            (
                IOobject
                (
                    "VReRes",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensioned<scalar>(word(),dimless,pTraits<scalar>::zero),
                fixedValueFvPatchField<scalar>::typeName
            );
            volScalarField VImRes
            (
                IOobject
                (
                    "VImRes",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensioned<scalar>(word(),dimless,pTraits<scalar>::zero),
                fixedValueFvPatchField<scalar>::typeName
            );

            VReRes.internalField().replace(0, VResIn.component(0));
            VImRes.internalField().replace(0, VResIn.component(1));

            VReRes.write();
            VImRes.write();
        }
    }

//     // Solve V only in conducting subregion
//     {
//         // Create submesh
//         fvMeshSubset cMeshSubSet
//         (
//             IOobject
//             (
//                 "conductor",
//                 mesh.time().timeName(),
//                 mesh,
//                 IOobject::NO_READ
//             ),
//             mesh
//         );
//
//         // Initialize large cell subset
//         cMeshSubSet.setLargeCellSubset(cCellsHashSet);
//
//         // Create sub-mesh
//         fvMesh& cMesh = cMeshSubSet.subMesh();
//
//         // References for maps
//         const labelList& cCellMap = cMeshSubSet.cellMap();
//
//         // Create fields
//         volScalarField cSigma
//         (
//             IOobject
//             (
//                 "sigma",
//                 cMesh.time().timeName(),
//                 cMesh,
//                 IOobject::NO_READ
//             ),
//             cMesh,
//             dimensionedScalar(word(), sigma.dimensions(), pTraits<scalar>::zero),
//             zeroGradientFvPatchField<scalar>::typeName
//         );
//
//         volScalarField cAlpha
//         (
//             IOobject
//             (
//                 "alpha",
//                 cMesh.time().timeName(),
//                 cMesh,
//                 IOobject::NO_READ
//             ),
//             cMesh,
//             dimensionedScalar(word(), alpha.dimensions(), pTraits<scalar>::zero),
//             zeroGradientFvPatchField<scalar>::typeName
//         );
//
//         volVectorField cARe
//         (
//             IOobject
//             (
//                 "ARe",
//                 cMesh.time().timeName(),
//                 cMesh,
//                 IOobject::NO_READ
//             ),
//             cMesh,
//             dimensioned<vector>(word(), ARe.dimensions(), pTraits<vector>::zero),
//             zeroGradientFvPatchField<vector>::typeName
//         );
//
//         volVectorField cAIm
//         (
//             IOobject
//             (
//                 "AIm",
//                 cMesh.time().timeName(),
//                 cMesh,
//                 IOobject::NO_READ
//             ),
//             cMesh,
//             dimensionedVector(word(), AIm.dimensions(), vector::zero),
//             zeroGradientFvPatchField<vector>::typeName
//         );
//
//         volScalarField cVRe
//         (
//             IOobject
//             (
//                 "VRe",
//                 cMesh.time().timeName(),
//                 cMesh,
//                 IOobject::NO_READ
//             ),
//             cMesh,
//             dimensionedScalar(word(), VRe.dimensions(), pTraits<scalar>::zero),
//             zeroGradientFvPatchField<scalar>::typeName
//         );
//
//         volScalarField cVIm
//         (
//             IOobject
//             (
//                 "VIm",
//                 cMesh.time().timeName(),
//                 cMesh,
//                 IOobject::NO_READ
//             ),
//             cMesh,
//             dimensionedScalar(word(), VIm.dimensions(), pTraits<scalar>::zero),
//             zeroGradientFvPatchField<scalar>::typeName
//         );
//
//         // Map internal values of volume fields
//         // and correct their boundary values
//         {
//             forAll (cMesh.C(), cellI)
//             {
//                 cSigma.internalField()[cellI] = sigma.internalField()[cCellMap[cellI]];
//                 cAlpha.internalField()[cellI] = alpha.internalField()[cCellMap[cellI]];
//
//                 cARe.internalField()[cellI] = ARe.internalField()[cCellMap[cellI]];
//                 cAIm.internalField()[cellI] = AIm.internalField()[cCellMap[cellI]];
//
//                 cVRe.internalField()[cellI] = VRe.internalField()[cCellMap[cellI]];
//                 cVIm.internalField()[cellI] = VIm.internalField()[cCellMap[cellI]];
//             }
//
//             cSigma.correctBoundaryConditions();
//             cAlpha.correctBoundaryConditions();
//
//             cARe.correctBoundaryConditions();
//             cAIm.correctBoundaryConditions();
//
//             cVRe.correctBoundaryConditions();
//             cVIm.correctBoundaryConditions();
//         }
//
//         // Solve V equation
//         {
//             // Assemble V-equations
//             fvScalarMatrix cVReEqn
//             (
//                 fvm::laplacian(cSigma, cVRe, "laplacian(sigmaV,VRe)")
//               - fvc::div(cAlpha*cAIm, "div(AIm)")
//             );
//             fvScalarMatrix cVImEqn
//             (
//                 fvm::laplacian(cSigma, cVIm, "laplacian(sigmaV,VIm)")
//               + fvc::div(cAlpha*cARe, "div(ARe)")
//             );
//
//             label refCell = cMesh.findCell(point(0.00, 0.00, 0.01));
//             scalar refValue = 0.0;
//
//             cVReEqn.setReference(refCell, refValue);
//             cVImEqn.setReference(refCell, refValue);
//
//             // Solve V
//             cVReEqn.solve();
//             cVImEqn.solve();
//
//             cVRe.correctBoundaryConditions();
//             cVIm.correctBoundaryConditions();
//         }
//
//         // Reverse map fields
//         {
//             forAll (cMesh.C(), cellI)
//             {
//                 VRe[cCellMap[cellI]] = cVRe[cellI];
//                 VIm[cCellMap[cellI]] = cVIm[cellI];
//             }
//         }
//     }
