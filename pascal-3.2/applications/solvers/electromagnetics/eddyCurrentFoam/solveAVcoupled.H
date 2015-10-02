
    {
        // Create scale factor for V-equation compared
        // to A-equation based on dimensional relation
        // dim[A] <=> dim[grad(V)]
        dimensionedScalar scaleV
        (
            "scaleV",
//             pow(average(mesh.V()), scalar(1.0/3.0))
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

        // Coupling terms for real/imagninary ARe-AIm
        // dependence of magnetic vector potential
        {
            AVEqn.insertEquationCoupling(0, 4,  alphaIn);
            AVEqn.insertEquationCoupling(4, 0, -alphaIn);
            AVEqn.insertEquationCoupling(1, 5,  alphaIn);
            AVEqn.insertEquationCoupling(5, 1, -alphaIn);
            AVEqn.insertEquationCoupling(2, 6,  alphaIn);
            AVEqn.insertEquationCoupling(6, 2, -alphaIn);
        }

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

            // Correct implicit gradient of V at
            // conductor interface
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

                    // Write
                    {
                        // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
//                         // Real part
//                         vector cReUf = diag(SnCmptMag[faceI] * VReInAReU[faceI]) & vector::one;
//                         vector cReLf = diag(SnCmptMag[faceI] * VReInAReL[faceI]) & vector::one;
//
//                         VReInAReD[own[faceI]] += cReUf;
//                         VReInAReD[nei[faceI]] += cReLf;
//                         VReInAReU[faceI] -= cReUf;
//                         VReInAReL[faceI] -= cReLf;
//
//                         // Imaginary part
//                         vector cImUf = diag(SnCmptMag[faceI] * VImInAImU[faceI]) & vector::one;
//                         vector cImLf = diag(SnCmptMag[faceI] * VImInAImL[faceI]) & vector::one;
//
//                         VImInAImD[own[faceI]] += cImUf;
//                         VImInAImD[nei[faceI]] += cImLf;
//                         VImInAImU[faceI] -= cImUf;
//                         VImInAImL[faceI] -= cImLf;
                    }
                }

//                 // Remove coeffs for non-conducting cells
//                 // TODO: In how far will this change upper/lower of
//                 // cells in contact with non-conducting region?
//                 VReInARe *= indicatorVc;
//                 VImInAIm *= indicatorVc;
            }

            // Insert implicit gradient terms
            AVEqn.insertBlockCoupling(0, 3, VReInARe, true);
            AVEqn.insertBlockCoupling(4, 7, VImInAIm, true);
        }

        // Coupling terms of magnetic vector potential A
        // in equations for electric scalar potential V
        {
            // Create BlockLduSystems
            BlockLduSystem<vector, scalar> AImInVRe(mesh);
            BlockLduSystem<vector, scalar> AReInVIm(mesh);

            // Set source fields of BlockLduSystems to zero
            AImInVRe.source() = pTraits<scalar>::zero;
            AReInVIm.source() = pTraits<scalar>::zero;

            // Grab parts of BlockLduSystems as linear
            CoeffField<vector>::linearTypeField& AImInVReD = AImInVRe.diag().asLinear();
            CoeffField<vector>::linearTypeField& AImInVReU = AImInVRe.upper().asLinear();
            CoeffField<vector>::linearTypeField& AImInVReL = AImInVRe.lower().asLinear();
            CoeffField<vector>::linearTypeField& AReInVImD = AReInVIm.diag().asLinear();
            CoeffField<vector>::linearTypeField& AReInVImU = AReInVIm.upper().asLinear();
            CoeffField<vector>::linearTypeField& AReInVImL = AReInVIm.lower().asLinear();

            // Implicit gradient constraints
            // at conductor interface as coefficients
            // of BlockLduSystem
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

                    // Write
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

                    // TODO: What about gradient correction?
                    // We have to modifiy the source term, too. But
                    // this depends on the type of correction. Hmm...
                }
            }

            // Implicit source terms due to variable sigma
            // in conducting region with implied condition: div(A)=0
            {
                AImInVReD -= scaleV.value() * indicatorVcIn * gradAlphaIn;
                AReInVImD += scaleV.value() * indicatorVcIn * gradAlphaIn;
            }

            AVEqn.insertBlockCoupling(3, 4, AImInVRe, false);
            AVEqn.insertBlockCoupling(7, 0, AReInVIm, false);
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

        // TODO: Remove after debug
        {
            vector8Field AVResIn = AVEqn.residual(AV.internalField());

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

            AReRes.internalField().replace(0, AVResIn.component(0));
            AReRes.internalField().replace(1, AVResIn.component(1));
            AReRes.internalField().replace(2, AVResIn.component(2));
            AImRes.internalField().replace(0, AVResIn.component(4));
            AImRes.internalField().replace(1, AVResIn.component(5));
            AImRes.internalField().replace(2, AVResIn.component(6));
            VReRes.internalField().replace(0, AVResIn.component(3));
            VImRes.internalField().replace(0, AVResIn.component(7));

            AReRes.write();
            AImRes.write();
            VReRes.write();
            VImRes.write();
        }
    }
