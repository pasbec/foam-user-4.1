
    // Solve A
    {
        // Coupling terms of electric scalar potential V
        // in equation for magnetic vector potential A
        {
            // Get relaxation factor
            scalar relax = AVmethodDict.lookupOrDefault<scalar>("relax",1.0);

            if (relax < (1.0-SMALL))
            {
                // Relaxation of gradient of V
                VReGrad.relax(relax);
                VImGrad.relax(relax);
            }

            // Explicit gradient terms of electric scalar potential V
            AReEqn -= sigma * VReGrad * indicatorVc;
            AImEqn -= sigma * VImGrad * indicatorVc;
        }

        // References to internal alpha field
        const scalarField& alphaIn = alpha.internalField();

        // Prepare block system for A
        fvBlockMatrix<vector6> AEqn(A);

        // Insert equations for A
        AEqn.insertEquation(0, AReEqn);
        AEqn.insertEquation(3, AImEqn);

        // Coupling terms for real/imagninary ARe-AIm
        // dependence of magnetic vector potential
        {
            AEqn.insertEquationCoupling(0, 3,  alphaIn);
            AEqn.insertEquationCoupling(3, 0, -alphaIn);
            AEqn.insertEquationCoupling(1, 4,  alphaIn);
            AEqn.insertEquationCoupling(4, 1, -alphaIn);
            AEqn.insertEquationCoupling(2, 5,  alphaIn);
            AEqn.insertEquationCoupling(5, 2, -alphaIn);
        }

        // Solve A
        AEqn.solve();

        // Retrieve solution
        AEqn.retrieveSolution(0, ARe.internalField());
        AEqn.retrieveSolution(3, AIm.internalField());

        ARe.correctBoundaryConditions();
        AIm.correctBoundaryConditions();

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
    }

    // Solve V
    {
        // Coupling terms of magnetic vector potential A
        // in equation for electric scalar potential V
        {
            // Explicit gradient constraints
            // at conductor interface
            {
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

                    // Write source which represent gradient
                    // of V depending on A
                    {
                        VReEqn.source()[own[faceI]] -= e * wP * (Sf[faceI] & AIm[own[faceI]]);
                        VReEqn.source()[own[faceI]] -= e * wN * (Sf[faceI] & AIm[nei[faceI]]);
                        VReEqn.source()[nei[faceI]] += e * wP * (Sf[faceI] & AIm[own[faceI]]);
                        VReEqn.source()[nei[faceI]] += e * wN * (Sf[faceI] & AIm[nei[faceI]]);

                        VImEqn.source()[own[faceI]] += e * wP * (Sf[faceI] & ARe[own[faceI]]);
                        VImEqn.source()[own[faceI]] += e * wN * (Sf[faceI] & ARe[nei[faceI]]);
                        VImEqn.source()[nei[faceI]] -= e * wP * (Sf[faceI] & ARe[own[faceI]]);
                        VImEqn.source()[nei[faceI]] -= e * wN * (Sf[faceI] & ARe[nei[faceI]]);
                    }

                    // TODO: What about gradient correction?
                    // We have to modifiy the source term, too. But
                    // this depends on the type of correction. Hmm...
                }
            }

            // Explicit source terms due to variable sigma
            // in conducting region with implied condition div(A)=0
            VReEqn -= (gradAlpha & AIm) * indicatorVc;
            VImEqn += (gradAlpha & ARe) * indicatorVc;
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

//         // TODO: Remove after debug
//         {
//             vector2Field VResIn = VEqn.residual(V.internalField());
//
//             volScalarField VReRes
//             (
//                 IOobject
//                 (
//                     "VReRes",
//                     mesh.time().timeName(),
//                     mesh,
//                     IOobject::NO_READ,
//                     IOobject::NO_WRITE
//                 ),
//                 mesh,
//                 dimensioned<scalar>(word(),dimless,pTraits<scalar>::zero),
//                 fixedValueFvPatchField<scalar>::typeName
//             );
//             volScalarField VImRes
//             (
//                 IOobject
//                 (
//                     "VImRes",
//                     mesh.time().timeName(),
//                     mesh,
//                     IOobject::NO_READ,
//                     IOobject::NO_WRITE
//                 ),
//                 mesh,
//                 dimensioned<scalar>(word(),dimless,pTraits<scalar>::zero),
//                 fixedValueFvPatchField<scalar>::typeName
//             );
//
//             VReRes.internalField().replace(0, VResIn.component(0));
//             VImRes.internalField().replace(0, VResIn.component(1));
//
//             VReRes.write();
//             VImRes.write();
//         }
    }
