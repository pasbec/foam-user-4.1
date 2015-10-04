
    {
        // Create scale factor for V-equation compared
        // to A-equation based on dimensional relation
        // dim[A] <=> dim[grad(V)]
        dimensionedScalar scaleV
        (
            "scaleV",
            pow(min(mesh.V()), scalar(1.0/3.0))
        );

        // Multiply V-equations with scale factor
        VReEqn *= scaleV;
        VImEqn *= scaleV;

        // Prepare block system
        fvBlockMatrix<vector8> AVEqn(AV);

        // Insert equations for A
        AVEqn.insertEquation(0, AReEqn);
        AVEqn.insertEquation(4, AImEqn);

        // Insert equations for V
        AVEqn.insertEquation(3, VReEqn);
        AVEqn.insertEquation(7, VImEqn);

        // References to internal fields
        const scalarField& sigmaIn = sigma.internalField();
        const scalarField& alphaIn = alpha.internalField();
        const vectorField& gradAlphaIn = gradAlpha.internalField();
        const scalarField& indicatorVcIn = indicatorVc.internalField();

        // Coupling terms of electric scalar potential V
        // in equation for magnetic vector potential A
        {
            // Implicit gradient terms of electric scalar potential V
            BlockLduSystem<vector, vector> VReInARe(fvm::grad(VRe)); VReInARe *= -sigmaIn;
            BlockLduSystem<vector, vector> VImInAIm(fvm::grad(VIm)); VImInAIm *= -sigmaIn;

            // Grab parts of BlockLduSystem as linear
            CoeffField<vector>::linearTypeField& VReInAReD = VReInARe.diag().asLinear();
            CoeffField<vector>::linearTypeField& VReInAReU = VReInARe.upper().asLinear();
            CoeffField<vector>::linearTypeField& VReInAReL = VReInARe.lower().asLinear();
            CoeffField<vector>::linearTypeField& VImInAImD = VImInAIm.diag().asLinear();
            CoeffField<vector>::linearTypeField& VImInAImU = VImInAIm.upper().asLinear();
            CoeffField<vector>::linearTypeField& VImInAImL = VImInAIm.lower().asLinear();

            // Force zero gradient of V (dV/dn=0)
            // at conductor interface as inner constraint.. We
            // need to do this in order set a new gradient
            // for V depnding on A, later!
            {
                // Mesh cell-face addressing
                const unallocLabelList& own = mesh.owner();
                const unallocLabelList& nei = mesh.neighbour();

                // Mesh geometry
                const vectorField& Sf = mesh.Sf();

                // Loop over all faces of conductor interface
                forAll (cnFaces, facei)
                {
                    label faceI = cnFaces[facei];

                    // Remove all corresponding coeffs from
                    // implicit face interpolation from owner
                    // and neighbour
                    {
                        VReInAReD[own[faceI]] += VReInAReU[faceI];
                        VReInAReD[nei[faceI]] += VReInAReL[faceI];
                        VReInAReU[faceI]      -= VReInAReU[faceI];
                        VReInAReL[faceI]      -= VReInAReL[faceI];

                        VImInAImD[own[faceI]] += VImInAImU[faceI];
                        VImInAImD[nei[faceI]] += VImInAImL[faceI];
                        VImInAImU[faceI]      -= VImInAImU[faceI];
                        VImInAImL[faceI]      -= VImInAImL[faceI];
                    }

                    // Write diag coeffs to set zero gradient
                    {
                        VReInAReD[own[faceI]] += Sf[faceI];
                        VReInAReD[nei[faceI]] -= Sf[faceI];

                        VImInAImD[own[faceI]] += Sf[faceI];
                        VImInAImD[nei[faceI]] -= Sf[faceI];
                    }
                }
            }

            // Insert implicit gradient terms
            AVEqn.insertBlockCoupling(0, 3, VReInARe);
            AVEqn.insertBlockCoupling(4, 7, VImInAIm);
        }

        // Coupling terms for real/imagninary parts of A
        {
            // Explicitely create empty BlockLduSystems
            BlockLduSystem<tensor, vector> AImInARe(mesh);
            BlockLduSystem<tensor, vector> AReInAIm(mesh);

            // Grab parts of BlockLduSystems as linear
            CoeffField<tensor>::linearTypeField& AImInAReD = AImInARe.diag().asLinear();
            CoeffField<tensor>::linearTypeField& AImInAReU = AImInARe.upper().asLinear();
            CoeffField<tensor>::linearTypeField& AImInAReL = AImInARe.lower().asLinear();
            CoeffField<tensor>::linearTypeField& AReInAImD = AReInAIm.diag().asLinear();
            CoeffField<tensor>::linearTypeField& AReInAImU = AReInAIm.upper().asLinear();
            CoeffField<tensor>::linearTypeField& AReInAImL = AReInAIm.lower().asLinear();

            // Point-implicit coupling from complex
            // harmonic time derivatives
            {
                // Mesh geometry
                const scalarField& V = mesh.V();

                AImInAReD += tensor(1,0,0,0,1,0,0,0,1) * alphaIn * V;
                AReInAImD -= tensor(1,0,0,0,1,0,0,0,1) * alphaIn * V;
            }

            AVEqn.insertBlockCoupling(0, 4, AImInARe);
            AVEqn.insertBlockCoupling(4, 0, AReInAIm);
        }

        // Coupling terms of magnetic vector potential A
        // in equations for electric scalar potential V
        {
            // Explicitely create empty BlockLduSystems
            BlockLduSystem<vector, scalar> AImInVRe(mesh);
            BlockLduSystem<vector, scalar> AReInVIm(mesh);

            // Grab parts of BlockLduSystems as linear
            CoeffField<vector>::linearTypeField& AImInVReD = AImInVRe.diag().asLinear();
            CoeffField<vector>::linearTypeField& AImInVReU = AImInVRe.upper().asLinear();
            CoeffField<vector>::linearTypeField& AImInVReL = AImInVRe.lower().asLinear();
            CoeffField<vector>::linearTypeField& AReInVImD = AReInVIm.diag().asLinear();
            CoeffField<vector>::linearTypeField& AReInVImU = AReInVIm.upper().asLinear();
            CoeffField<vector>::linearTypeField& AReInVImL = AReInVIm.lower().asLinear();

            // Implicit gradient constraints at conductor
            // interface as coefficients of BlockLduSystem
            {
                // Mesh cell-face addressing
                const unallocLabelList& own = mesh.owner();
                const unallocLabelList& nei = mesh.neighbour();

                // Mesh geometry
                const vectorField& Cf = mesh.Cf();
                const vectorField& C = mesh.C();
                const vectorField& Sf = mesh.Sf();
                const scalarField& magSf = mesh.magSf();
                const vectorField Sn = Sf/magSf;

                // Loop over all faces of conductor interface
                forAll (cnFaces, facei)
                {
                    label faceI = cnFaces[facei];

                    // Calculate face centre to cell centre distances Pf and fN
                    scalar Pf = mag(Sn[faceI] & (Cf[faceI] - C[own[faceI]]));
                    scalar fN = mag(Sn[faceI] & (C[nei[faceI]] - Cf[faceI]));

                    // Calculate weights
                    scalar wP = fN / (Pf+fN);
                    scalar wN = Pf / (Pf+fN);

                    // Calculate gradient coeffs
                    scalar cP = scaleV.value() * omega.value() * sSigmaV[faceI] * magSf[faceI] * wP;
                    scalar cN = scaleV.value() * omega.value() * sSigmaV[faceI] * magSf[faceI] * wN;

                    // Write coeffs which represent gradient
                    // of V depending on A
                    {
                        AImInVReD[own[faceI]] += cP * Sn[faceI];
                        AImInVReU[faceI]      += cN * Sn[faceI];
                        AImInVReL[faceI]      -= cP * Sn[faceI];
                        AImInVReD[nei[faceI]] -= cN * Sn[faceI];

                        AReInVImD[own[faceI]] -= cP * Sn[faceI];
                        AReInVImU[faceI]      -= cN * Sn[faceI];
                        AReInVImL[faceI]      += cP * Sn[faceI];
                        AReInVImD[nei[faceI]] += cN * Sn[faceI];
                    }
                }
            }

            // Implicit source terms due to variable sigma
            // in conducting region with implied condition div(A)=0
            {
                AImInVReD -= scaleV.value() * indicatorVcIn * gradAlphaIn;
                AReInVImD += scaleV.value() * indicatorVcIn * gradAlphaIn;
            }

            AVEqn.insertBlockCoupling(3, 4, AImInVRe);
            AVEqn.insertBlockCoupling(7, 0, AReInVIm);
        }

        // Solve A
        AVEqn.solve();

        // Retrieve solution for A
        AVEqn.retrieveSolution(0, ARe.internalField());
        AVEqn.retrieveSolution(4, AIm.internalField());

        ARe.correctBoundaryConditions();
        AIm.correctBoundaryConditions();

        // Retrieve solution for V
        AVEqn.retrieveSolution(3, VRe.internalField());
        AVEqn.retrieveSolution(7, VIm.internalField());

        VRe.correctBoundaryConditions();
        VIm.correctBoundaryConditions();

//         // TODO: Remove after debug
//         {
//             vector8Field AVResIn = AVEqn.residual(AV.internalField());
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
//             AReRes.internalField().replace(0, AVResIn.component(0));
//             AReRes.internalField().replace(1, AVResIn.component(1));
//             AReRes.internalField().replace(2, AVResIn.component(2));
//             AImRes.internalField().replace(0, AVResIn.component(4));
//             AImRes.internalField().replace(1, AVResIn.component(5));
//             AImRes.internalField().replace(2, AVResIn.component(6));
//             VReRes.internalField().replace(0, AVResIn.component(3));
//             VImRes.internalField().replace(0, AVResIn.component(7));
//
//             AReRes.write();
//             AImRes.write();
//             VReRes.write();
//             VImRes.write();
//         }
    }
