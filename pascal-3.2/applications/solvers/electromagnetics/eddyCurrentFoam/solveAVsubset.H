
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

        // Map volume fields
        forAll (cMesh.C(), cellI)
        {
            cARe.internalField()[cellI] = ARe.internalField()[cMeshCellMap[cellI]];
            cAIm.internalField()[cellI] = AIm.internalField()[cMeshCellMap[cellI]];

            cVRe.internalField()[cellI] = VRe.internalField()[cMeshCellMap[cellI]];
            cVIm.internalField()[cellI] = VIm.internalField()[cMeshCellMap[cellI]];

            cSigma.internalField()[cellI] = sigma.internalField()[cMeshCellMap[cellI]];
        }
        forAll(cMesh.boundaryMesh(), patchI)
        {
            polyPatch patch = cMesh.boundaryMesh()[patchI];

            // Map old internal faces
            if (patch.name() == "oldInternalFaces")
            {
                forAll(patch, facei)
                {
                    label faceI = patch.start() + facei;

                    cARe.boundaryField()[patchI][facei] == sARe[cMeshFaceMap[faceI]];
                    cAIm.boundaryField()[patchI][facei] == sAIm[cMeshFaceMap[faceI]];
                }
            }

            // TODO: What about other patches?
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
                cVReP.gradient() =  omega.value()
                   * (cAIm.boundaryField()[patchI] & cVReP.patch().nf());
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
                cVImP.gradient() = -omega.value()
                   * (cARe.boundaryField()[patchI] & cVImP.patch().nf());
            }
        }

        // Interpolate sigma to face centers
        surfaceScalarField csSigma
        (
            "sSigma",
            fvc::interpolate(cSigma, "interpolate(sigma)")
        );

        // Prepare V-equations
        fvScalarMatrix cVReEqn(fvm::laplacian(csSigma, cVRe, "laplacian(sigma,VRe)"));
        fvScalarMatrix cVImEqn(fvm::laplacian(csSigma, cVIm, "laplacian(sigma,VIm)"));

        volVectorField cSigmaGrad = fvc::grad(cSigma);

        // Explicit source terms due to variable sigma
        // in conducting region with implied condition div(A)=0
        volScalarField cVReEqnSource =   omega * (cSigmaGrad & cAIm);
        volScalarField cVImEqnSource = - omega * (cSigmaGrad & cARe);

//         // TEST
//         cVReEqnSource +=   omega * cSigma * fvc::div(cAIm);
//         cVImEqnSource += - omega * cSigma * fvc::div(cARe);

        cVReEqn -= cVReEqnSource;
        cVImEqn -= cVImEqnSource;

        // Check Poisson compatibility
        {
            // Sums
            scalar cVReEqnSumSource = 0.0;
            scalar cVImEqnSumSource = 0.0;

            // Mesh geometry
            const scalarField& V = cMesh.V();

            forAll (cMesh.C(), cellI)
            {
                cVReEqnSumSource += cVReEqnSource[cellI] * V[cellI];
                cVImEqnSumSource += cVImEqnSource[cellI] * V[cellI];
            }

            scalar cVReEqnSumFlux = 0.0;
            scalar cVImEqnSumFlux = 0.0;

            forAll(cMesh.boundaryMesh(), patchI)
            {
                const polyPatch& patch = cMesh.boundaryMesh()[patchI];

                vectorField Nf = cMesh.Sf().boundaryField()[patchI] / cMesh.magSf().boundaryField()[patchI];

                // Old internal faces
                if (patch.name() == "oldInternalFaces")
                {
                    forAll(patch, facei)
                    {
                        label faceI = patch.start() + facei;

                        cVReEqnSumFlux += omega.value() * (sAIm[cMeshFaceMap[faceI]] & Nf[cMeshFaceMap[faceI]]);
                        cVImEqnSumFlux += omega.value() * (sARe[cMeshFaceMap[faceI]] & Nf[cMeshFaceMap[faceI]]);
                    }
                }

                // TODO: What about other patches?
            }

            Info << "VRe Poisson equation compatibility check (Source, Flux):"
                << "(" << cVReEqnSumSource << ", " << cVReEqnSumFlux << ")" << endl;
            Info << "VIm Poisson equation compatibility check (Source, Flux):"
                << "(" << cVImEqnSumSource << ", " << cVImEqnSumFlux << ")" << endl;
        }

        // Reference settings
        {
            // Reference values in conducting region
            scalar cRefValueVRe = 0.0;
            scalar cRefValueVIm = 0.0;

            // Reference cells in conducting region
            // FIXME: Refpoint problems!
            cVReEqn.setReference(csRefCell, cRefValueVRe);
            cVImEqn.setReference(csRefCell, cRefValueVIm);
        }

        // Solve V-equations
        cVReEqn.solve();
        cVImEqn.solve();

        cVRe.correctBoundaryConditions();
        cVIm.correctBoundaryConditions();

        // FIXME: Wrong gradient values at boundary
        volVectorField cVReGrad
        (
            IOobject
            (
                "VRe",
                cMesh.time().timeName(),
                cMesh,
                IOobject::NO_READ
            ),
            fvc::grad(cVRe)
        );

        // FIXME: Wrong gradient values at boundary
        volVectorField cVImGrad
        (
            IOobject
            (
                "VIm",
                cMesh.time().timeName(),
                cMesh,
                IOobject::NO_READ
            ),
            fvc::grad(cVIm)
        );

        // Reverse map V
        forAll (cMesh.C(), cellI)
        {
            VRe[cMeshCellMap[cellI]] = cVRe[cellI];
            VIm[cMeshCellMap[cellI]] = cVIm[cellI];

            VReGrad[cMeshCellMap[cellI]] = cVReGrad[cellI];
            VImGrad[cMeshCellMap[cellI]] = cVImGrad[cellI];
        }
    }

    // TEST
    // FIXME: Wrong gradient values at boundary (see above)
    // Update gradient of V
    {
        VReGrad = fvc::grad(VRe);
        VImGrad = fvc::grad(VIm);
    }
