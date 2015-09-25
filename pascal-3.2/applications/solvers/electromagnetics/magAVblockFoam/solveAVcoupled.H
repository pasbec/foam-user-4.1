
    {
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

        // Prepare block system
        fvBlockMatrix<vector8> AVEqn(AV);

        // Insert equations for A
        AVEqn.insertEquation(0, ARexEqn);
        AVEqn.insertEquation(1, AReyEqn);
        AVEqn.insertEquation(2, ARezEqn);
        AVEqn.insertEquation(4, AImxEqn);
        AVEqn.insertEquation(5, AImyEqn);
        AVEqn.insertEquation(6, AImzEqn);

        // Insert equations for V
        AVEqn.insertEquation(3, VReEqn);
        AVEqn.insertEquation(7, VImEqn);

        // References to internal fields
        const scalarField& sigmaIn = sigma.internalField();
        const scalarField& alphaIn = alpha.internalField();

        // Coupling terms for real/imagninary ARe-AIm coupling
        // of magnetic vector potential
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

            // Implicit gradient operators
            BlockLduSystem<vector, vector> VReInARe(fvm::grad(VRe));
            BlockLduSystem<vector, vector> VImInAIm(fvm::grad(VIm));

            // Term coefficients
            VReInARe *= -sigmaIn;
            VImInAIm *= -sigmaIn;

            // Modify fvm::grad(VRe)/fvm::grad(VIm) to enforce zero normal
            // gradient for V (d(V)/dn=0) at conductor interface as inner
            // constraint
            // TODO: Gradient correction for div(j)=0???
            {
                // Mesh addressing
                const unallocLabelList& own = mesh.owner();
                const unallocLabelList& nei = mesh.neighbour();

                // Mesh geometry
                const vectorField& Sf = mesh.Sf();
                const scalarField& magSf = mesh.magSf();

                // Derived mesh geometry
                const vectorField Sn = Sf/magSf;
                const vectorField SnCmptMag = cmptMag(Sn);

                // Reference to parts of blockLduSystem
                CoeffField<vector>::linearTypeField& VReInAReD = VReInARe.diag().asLinear();
                CoeffField<vector>::linearTypeField& VReInAReU = VReInARe.upper().asLinear();
                CoeffField<vector>::linearTypeField& VReInAReL = VReInARe.lower().asLinear();

                CoeffField<vector>::linearTypeField& VImInAImD = VImInAIm.diag().asLinear();
                CoeffField<vector>::linearTypeField& VImInAImU = VImInAIm.upper().asLinear();
                CoeffField<vector>::linearTypeField& VImInAImL = VImInAIm.lower().asLinear();

                // Loop over all faces of conductor interface
                forAll (cnFaces, facei)
                {
                    label faceI = cnFaces[facei];

                    // Eliminate upper/lower coeffs based on constraint
                    // NOTE: This will keep SYMMETRIC matrices and it will
                    //       also sligthly improve their diagonal dominance
                    {
                        // Real part
                        vector cReUf = diag(SnCmptMag[faceI] * VReInAReU[faceI]) & vector::one;
                        vector cReLf = diag(SnCmptMag[faceI] * VReInAReL[faceI]) & vector::one;

                        VReInAReD[own[faceI]] += cReUf;
                        VReInAReD[nei[faceI]] += cReLf;
                        VReInAReU[faceI] -= cReUf;
                        VReInAReL[faceI] -= cReLf;

                        // Imaginary part
                        vector cImUf = diag(SnCmptMag[faceI] * VImInAImU[faceI]) & vector::one;
                        vector cImLf = diag(SnCmptMag[faceI] * VImInAImL[faceI]) & vector::one;

                        VImInAImD[own[faceI]] += cImUf;
                        VImInAImD[nei[faceI]] += cImLf;
                        VImInAImU[faceI] -= cImUf;
                        VImInAImL[faceI] -= cImLf;
                    }
                }

                // Remove coeffs for non-conducting cells
                // TODO: In how far will this change upper/lower of
                // cells in contact with non-conducting region?
                VReInARe *= sigmaIc;
                VImInAIm *= sigmaIc;
            }

            // Insert coupling of V in A
            AVEqn.insertBlockCoupling(0, 3, VReInARe, true);
            AVEqn.insertBlockCoupling(4, 7, VImInAIm, true);
        }

        // Coupling terms of magnetic vector potential A
        // in equation for electric scalar potential V
        {
            // Implicit divergence operators
            BlockLduSystem<vector, scalar> AImInVRe(fvm::UDiv(sSigmaV, AIm, "div(AIm)"));
            BlockLduSystem<vector, scalar> AReInVIm(fvm::UDiv(sSigmaV, ARe, "div(ARe)"));

            AImInVRe *= -scaleV.value() * omega.value();
            AReInVIm *=  scaleV.value() * omega.value();

            // Modify fvm::UDiv(ARe)/fvm::UDiv(AIm) to enforce zero
            // normal components (ARe*n=0) and (AIm*n=0) at the
            // conductor interface as an inner constraint
            {
                // Mesh addressing
                const unallocLabelList& own = mesh.owner();
                const unallocLabelList& nei = mesh.neighbour();

                // Mesh geometry
                const vectorField& Sf = mesh.Sf();
                const scalarField& magSf = mesh.magSf();

                // Derived mesh geometry
                const vectorField Sn = Sf/magSf;
                const vectorField SnCmptMag = cmptMag(Sn);

                // Reference to parts of blockLduSystem
                CoeffField<vector>::linearTypeField& AImInVReD = AImInVRe.diag().asLinear();
                CoeffField<vector>::linearTypeField& AImInVReU = AImInVRe.upper().asLinear();
                CoeffField<vector>::linearTypeField& AImInVReL = AImInVRe.lower().asLinear();

                CoeffField<vector>::linearTypeField& AReInVImD = AReInVIm.diag().asLinear();
                CoeffField<vector>::linearTypeField& AReInVImU = AReInVIm.upper().asLinear();
                CoeffField<vector>::linearTypeField& AReInVImL = AReInVIm.lower().asLinear();

                // Loop over all faces of conductor interface
                forAll (cnFaces, facei)
                {
                    label faceI = cnFaces[facei];

                    // Eliminate upper/lower coeffs based on constraint
                    // NOTE: This will keep SYMMETRIC matrices and it will
                    //       also sligthly improve their diagonal dominance
                    {
                        // Real part
                        vector cReUf = diag(SnCmptMag[faceI] * AImInVReU[faceI]) & vector::one;
                        vector cReLf = diag(SnCmptMag[faceI] * AImInVReL[faceI]) & vector::one;

                        AImInVReD[own[faceI]] += cReLf;
                        AImInVReD[nei[faceI]] += cReUf;
                        AImInVReU[faceI] -= cReUf;
                        AImInVReL[faceI] -= cReLf;

                        // Imaginary part
                        vector cImUf = diag(SnCmptMag[faceI] * AReInVImU[faceI]) & vector::one;
                        vector cImLf = diag(SnCmptMag[faceI] * AReInVImL[faceI]) & vector::one;

                        AReInVImD[own[faceI]] += cImLf;
                        AReInVImD[nei[faceI]] += cImUf;
                        AReInVImU[faceI] -= cImUf;
                        AReInVImL[faceI] -= cImLf;
                    }
                }

                // Remove coeffs for non-conducting cells
                // TODO: In how far will this change upper/lower of
                // cells in contact with non-conducting region?
                AImInVRe *= sigmaIc;
                AReInVIm *= sigmaIc;
            }

            // Insert coupling of A in V
            AVEqn.insertBlockCoupling(3, 4, AImInVRe, false);
            AVEqn.insertBlockCoupling(7, 0, AReInVIm, false);
        }

        // Solve A
        AVEqn.solve();

        // Retrieve solution for A
        AVEqn.retrieveSolution(0, ARex.internalField());
        AVEqn.retrieveSolution(1, ARey.internalField());
        AVEqn.retrieveSolution(2, ARez.internalField());
        AVEqn.retrieveSolution(4, AImx.internalField());
        AVEqn.retrieveSolution(5, AImy.internalField());
        AVEqn.retrieveSolution(6, AImz.internalField());

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

        // TODO: This needs to be the same as the
        // component boundary conditions otherwise
        // it gives non-sense. This means we currently
        // cannot use directionMixed bc!!!
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
