
    // Solve V
    {
        // Explicit gradient constraints
        // at conductor interface
        {
            // Grab source of V-equations
            scalarField& VReEqnS = VReEqn.source();
            scalarField& VImEqnS = VImEqn.source();

            // Mesh cell-face addressing
            const unallocLabelList& own = mesh.owner();
            const unallocLabelList& nei = mesh.neighbour();

            // Mesh geometry
            const vectorField& Sf = mesh.Sf();
            const scalarField& w = mesh.weights();

            // Loop over all faces of conductor interface
            forAll (cnFaces, facei)
            {
                label faceI = cnFaces[facei];

                // Calculate weight coeffs
                scalar wP = omega.value() * w[faceI];
                scalar wN = omega.value() - wP;

                // Calculate equation coeffs
                scalar e = sSigmaV[faceI];

                // Write source which represents gradient
                // of V depending on A
                {
                    VReEqnS[own[faceI]] -= e * wP * (Sf[faceI] & AIm[own[faceI]]);
                    VReEqnS[own[faceI]] -= e * wN * (Sf[faceI] & AIm[nei[faceI]]);
                    VReEqnS[nei[faceI]] += e * wP * (Sf[faceI] & AIm[own[faceI]]);
                    VReEqnS[nei[faceI]] += e * wN * (Sf[faceI] & AIm[nei[faceI]]);

                    VImEqnS[own[faceI]] += e * wP * (Sf[faceI] & ARe[own[faceI]]);
                    VImEqnS[own[faceI]] += e * wN * (Sf[faceI] & ARe[nei[faceI]]);
                    VImEqnS[nei[faceI]] -= e * wP * (Sf[faceI] & ARe[own[faceI]]);
                    VImEqnS[nei[faceI]] -= e * wN * (Sf[faceI] & ARe[nei[faceI]]);
                }

                // TODO: What about gradient correction?
                // We have to modifiy the source term, too. But
                // this depends on the type of correction. Hmm...
            }
        }

        // Explicit source terms due to variable sigma
        // in conducting region with implied condition div(A)=0
        volScalarField VReEqnSource("VReEqnSource",   gradAlpha & AIm);
        volScalarField VImEqnSource("VImEqnSource", - gradAlpha & ARe);

//         // TEST
//         VReEqnSource +=   alpha * fvc::div(AIm);
//         VImEqnSource += - alpha * fvc::div(ARe);

//         // FIXME: Refpoint problems!
//         // Check/Correct Poisson compatibility
//         {
//             // Source sums
//             scalar VReEqnSumSource = 0.0;
//             scalar VImEqnSumSource = 0.0;
//
//             // Flux sums
//             scalar VReEqnSumFlux = 0.0;
//             scalar VImEqnSumFlux = 0.0;
//
//             // Correct source integrals
//             // TEST
//             {
//                 // Global conducting cell size
//                 label cCellsSizeGlobal = returnReduce(cCells.size(), sumOp<label>());
//
//                 // Global source integrals
//                 VReEqnSumSource = gSum(VReEqnSource.internalField() * mesh.V());
//                 VImEqnSumSource = gSum(VImEqnSource.internalField() * mesh.V());
//
//                 // Force global source integral to zero
//                 forAll (cCells, celli)
//                 {
//                     label cellI = cCells[celli];
//
//                     VReEqnSource[cellI] -= VReEqnSumSource / cCellsSizeGlobal / mesh.V()[cellI];
//                     VImEqnSource[cellI] -= VImEqnSumSource / cCellsSizeGlobal / mesh.V()[cellI];
//                 }
//             }
//
//             // Source integrals
//             VReEqnSumSource = gSum(VReEqnSource.internalField()*mesh.V());
//             VImEqnSumSource = gSum(VImEqnSource.internalField()*mesh.V());
//
//             // Boundary integrals
//             {
//                 surfaceScalarField SfARe("SfARe", linearInterpolate(ARe) & mesh.Sf());
//                 surfaceScalarField SfAIm("SfAIm", linearInterpolate(AIm) & mesh.Sf());
//
//                 forAll (cnPNfaces, facei)
//                 {
//                     label faceI = cnPNfaces[facei];
//
//                     VReEqnSumFlux += omega.value() * SfAIm[faceI];
//                     VImEqnSumFlux -= omega.value() * SfARe[faceI];
//                 }
//
//                 forAll (cnNPfaces, facei)
//                 {
//                     label faceI = cnNPfaces[facei];
//
//                     VReEqnSumFlux -= omega.value() * SfAIm[faceI];
//                     VImEqnSumFlux += omega.value() * SfARe[faceI];
//                 }
//
//                 reduce(VReEqnSumFlux, sumOp<scalar>());
//                 reduce(VImEqnSumFlux, sumOp<scalar>());
//             }
//
//             Info << "VRe Poisson equation compatibility check (Source, Flux):"
//                 << "(" << VReEqnSumSource << ", " << VReEqnSumFlux << ")" << endl;
//             Info << "VIm Poisson equation compatibility check (Source, Flux):"
//                 << "(" << VImEqnSumSource << ", " << VImEqnSumFlux << ")" << endl;
//         }

        // Add source terms
        VReEqn -= VReEqnSource;
        VImEqn -= VImEqnSource;

        // Reference points
        // FIXME: Refpoint problems!
        {
            // Grab parts of V-equations
            scalarField& VReEqnD = VReEqn.diag();
            scalarField& VReEqnS = VReEqn.source();
            scalarField& VImEqnD = VImEqn.diag();
            scalarField& VImEqnS = VImEqn.source();

            // Reference values in conducting region
            scalar cRefValueVRe = 0.0;
            scalar cRefValueVIm = 0.0;

            // Force reference cells in conducting region
            VReEqnD[cRefCell] += VReEqnD[cRefCell];
            VImEqnD[cRefCell] += VImEqnD[cRefCell];
            VReEqnS[cRefCell] += VReEqnD[cRefCell]*cRefValueVRe;
            VImEqnS[cRefCell] += VImEqnD[cRefCell]*cRefValueVIm;

            // Reference values in non-conducting region
            // TODO: Which values to use here best?
            scalar nRefValueVRe = 0.0;
            scalar nRefValueVIm = 0.0;

            // Set reference cells in non-conducting region if needed
            VReEqn.setReference(nRefCell, nRefValueVRe);
            VImEqn.setReference(nRefCell, nRefValueVIm);
        }

//         // Pseudo time-stepping based on Euler-explicit
//         // TEST
//         {
//             // Grab parts of V-equations
//             scalarField& VReEqnD = VReEqn.diag();
//             scalarField& VReEqnS = VReEqn.source();
//             scalarField& VImEqnD = VImEqn.diag();
//             scalarField& VImEqnS = VImEqn.source();
//
// //             // Fictive time step
// //             // NOTE: Minus sign as laplacian is in matrix with plus sign
// //             scalarField V = mesh.V();
// //             scalar pseudoDeltaT = 1.0;
// //             VReEqnD -= V/pseudoDeltaT;
// //             VImEqnD -= V/pseudoDeltaT;
// //             VReEqnS -= V/pseudoDeltaT * VRe.internalField();
// //             VImEqnS -= V/pseudoDeltaT * VIm.internalField();
//
//             // Arbitrary damp factor based on diagonal value
//             scalar damp = 1.0e-04;
//             VReEqnD += damp * VReEqnD;
//             VImEqnD += damp * VImEqnD;
//             VReEqnS += damp * VReEqnD * VRe.internalField();
//             VImEqnS += damp * VImEqnD * VIm.internalField();
//         }

        // Restore solution
        VRe += VReAverage;
        VIm += VImAverage;

//         // Solve V-equations
//         // FIXME: Why does GAMG or amgSolver not work well?
//         VReEqn.solve();
//         VImEqn.solve();
//
//         VRe.correctBoundaryConditions();
//         VIm.correctBoundaryConditions();

        // Create complex vector2 field
        volVector2Field V
        (
            IOobject
            (
                "V",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector2(word(),dimless,vector2::zero)
        );

        // Init complex vector2 field
        V.replace(0, VRe);
        V.replace(1, VIm);

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

        // Calculate average
        VReAverage = average(VRe);
        VImAverage = average(VIm);

        // Normalize solution
        VRe -= VReAverage;
        VIm -= VImAverage;

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

        // Update gradient of V
        VReGrad = fvc::grad(VRe);
        VImGrad = fvc::grad(VIm);
    }
