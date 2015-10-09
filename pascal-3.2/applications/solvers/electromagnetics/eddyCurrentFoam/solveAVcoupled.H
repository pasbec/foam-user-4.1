
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
            BlockLduSystem<vector, vector> VReInARe(fvm::grad(VRe));
            BlockLduSystem<vector, vector> VImInAIm(fvm::grad(VIm));

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
            // for V depnding on A in a second step.
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

                    // Write zero gradient coeffs
                    {
                        VReInAReD[own[faceI]] += Sf[faceI];
                        VReInAReD[nei[faceI]] -= Sf[faceI];

                        VImInAImD[own[faceI]] += Sf[faceI];
                        VImInAImD[nei[faceI]] -= Sf[faceI];
                    }
                }
            }

            // Multiply with equation coeffs
            VReInARe *= - sigmaIn * indicatorVcIn;
            VImInAIm *= - sigmaIn * indicatorVcIn;

            // Insert implicit gradient terms
            AVEqn.insertBlockCoupling(0, 3, VReInARe);
            AVEqn.insertBlockCoupling(4, 7, VImInAIm);
        }

        // Coupling terms for real/imagninary parts of A
        {
            // Explicitly create empty BlockLduSystems
            BlockLduSystem<tensor, vector> AImInARe(mesh);
            BlockLduSystem<tensor, vector> AReInAIm(mesh);

            // Grab parts of BlockLduSystems as linear
            CoeffField<tensor>::linearTypeField& AImInAReD = AImInARe.diag().asLinear();
            CoeffField<tensor>::linearTypeField& AImInAReU = AImInARe.upper().asLinear();
            CoeffField<tensor>::linearTypeField& AImInAReL = AImInARe.lower().asLinear();
            CoeffField<tensor>::linearTypeField& AReInAImD = AReInAIm.diag().asLinear();
            CoeffField<tensor>::linearTypeField& AReInAImU = AReInAIm.upper().asLinear();
            CoeffField<tensor>::linearTypeField& AReInAImL = AReInAIm.lower().asLinear();

            // Correct gradient
            {
                // Mesh cell-face addressing
                const unallocLabelList& own = mesh.owner();
                const unallocLabelList& nei = mesh.neighbour();

                // Mesh geometry
                const vectorField& Sf = mesh.Sf();
                const scalarField& magSf = mesh.magSf();
                const scalarField& w = mesh.weights();
                const scalarField& d = mesh.deltaCoeffs();

                // Derived mesh geometry
                const vectorField Sn = Sf/magSf;
                tensorField SfSn = Sf * Sn;

                // Set correct gradient
                forAll (cnFaces, facei)
                {
                    label faceI = cnFaces[facei];

                    // Calculate face centre to cell centre distances Pf and fN
                    scalar fN = w[faceI] / d[faceI];
                    scalar Pf = 1.0 / d[faceI] - fN;

                    // Calculate weight coeffs
                    scalar wP = omega.value() * w[faceI];
                    scalar wN = omega.value() - wP;

                    // Write coeffs which represent gradient part in V
                    // from other complex part of A
                    {
                        AImInAReD[own[faceI]] += Pf * wP * SfSn[faceI];
                        AImInAReU[faceI]      += Pf * wN * SfSn[faceI];
                        AImInAReL[faceI]      += fN * wP * SfSn[faceI];
                        AImInAReD[nei[faceI]] += fN * wN * SfSn[faceI];

                        AReInAImD[own[faceI]] -= Pf * wP * SfSn[faceI];
                        AReInAImU[faceI]      -= Pf * wN * SfSn[faceI];
                        AReInAImL[faceI]      -= fN * wP * SfSn[faceI];
                        AReInAImD[nei[faceI]] -= fN * wN * SfSn[faceI];
                    }
                }

                // Multiply with equation coeffs
                AImInARe *= - sigmaIn * indicatorVcIn;
                AReInAIm *= - sigmaIn * indicatorVcIn;
            }

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
            // Explicitly create empty BlockLduSystems
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
                const vectorField& Sf = mesh.Sf();
                const scalarField& w = mesh.weights();

                // Loop over all faces of conductor interface
                forAll (cnFaces, facei)
                {
                    label faceI = cnFaces[facei];

                    // Calculate gradient coeffs
                    scalar wP = omega.value() * w[faceI];
                    scalar wN = omega.value() - wP;

                    // Calculate equation coeffs
                    scalar e = scaleV.value() * sSigmaV[faceI];

                    // Write coeffs which represent gradient
                    // of V depending on A
                    {
                        AImInVReD[own[faceI]] += e * wP * Sf[faceI];
                        AImInVReU[faceI]      += e * wN * Sf[faceI];
                        AImInVReL[faceI]      -= e * wP * Sf[faceI];
                        AImInVReD[nei[faceI]] -= e * wN * Sf[faceI];

                        AReInVImD[own[faceI]] -= e * wP * Sf[faceI];
                        AReInVImU[faceI]      -= e * wN * Sf[faceI];
                        AReInVImL[faceI]      += e * wP * Sf[faceI];
                        AReInVImD[nei[faceI]] += e * wN * Sf[faceI];
                    }

                    // TODO: What about gradient correction?
                    // We have to modifiy the source term, too. But
                    // this depends on the type of correction. Hmm...
                }
            }

            // Implicit source terms due to variable sigma
            // in conducting region with implied condition div(A)=0
            {
                // Mesh geometry
                const scalarField& V = mesh.V();

                AImInVReD -= scaleV.value() * gradAlphaIn * V * indicatorVcIn;
                AReInVImD += scaleV.value() * gradAlphaIn * V * indicatorVcIn;
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

    // Update gradient of V
    {
        VReGrad = fvc::grad(VRe);
        VImGrad = fvc::grad(VIm);
    }
