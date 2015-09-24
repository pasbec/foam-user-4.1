
    // Solve Div(A) in A
    {
        volScalarField AReDiv
        (
            IOobject
            (
                "AReDiv",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::div(ARe),
            fixedValueFvPatchField<scalar>::typeName
        );

        volScalarField AImDiv
        (
            IOobject
            (
                "AImDiv",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::div(AIm),
            fixedValueFvPatchField<scalar>::typeName
        );

        {
            fvScalarMatrix AReDivEqn
            (
                rMu0 * fvm::Sp(scalar(1.0),AReDiv)
//            - rMu0 * fvc::div(ARe)
            );

            // Currently only work-around
            // since diagonal insertion
            // in block-matrix does not
            // work correctly
            AReDivEqn.upper() = 0.0;

            fvScalarMatrix AImDivEqn
            (
                rMu0 * fvm::Sp(scalar(1.0),AImDiv)
//            - rMu0 * fvc::div(AIm)
            );

            // Currently only work-around
            // since diagonal insertion
            // in block-matrix does not
            // work correctly
            AImDivEqn.upper() = 0.0;

            // Insert equations for divA
            AEqn.insertEquation(6, AReDivEqn);
            AEqn.insertEquation(7, AImDivEqn);

            // Assemble coupling of divA in A
            BlockLduSystem<vector, vector> AReDivInARe(fvm::grad(AReDiv)); AReDivInARe *= -rMu0.value();
            BlockLduSystem<vector, vector> AImDivInAIm(fvm::grad(AImDiv)); AImDivInAIm *= -rMu0.value();

            // Modify fvm::grad() for A*n=0 at conductor interface
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
                CoeffField<vector>::linearTypeField& AReDivInAReD = AReDivInARe.diag().asLinear();
                CoeffField<vector>::linearTypeField& AReDivInAReU = AReDivInARe.upper().asLinear();
                CoeffField<vector>::linearTypeField& AReDivInAReL = AReDivInARe.lower().asLinear();

                CoeffField<vector>::linearTypeField& AImDivInAImD = AImDivInAIm.diag().asLinear();
                CoeffField<vector>::linearTypeField& AImDivInAImU = AImDivInAIm.upper().asLinear();
                CoeffField<vector>::linearTypeField& AImDivInAImL = AImDivInAIm.lower().asLinear();

                // Loop over all faces of conductor interface
                forAll (cnFaces, facei)
                {
                    label faceI = cnFaces[facei];

                    // Eliminate all coeffs
                    // TODO: Optimize
                    {
                        // Real part
                        vector cReUf = diag(SnCmptMag[faceI] * AReDivInAReU[faceI]) & vector::one;
                        vector cReLf = diag(SnCmptMag[faceI] * AReDivInAReL[faceI]) & vector::one;

                        AReDivInAReD[own[faceI]] += cReUf;
                        AReDivInAReD[nei[faceI]] += cReLf;
                        AReDivInAReU[faceI] -= cReUf;
                        AReDivInAReL[faceI] -= cReLf;

                        // Imaginary part
                        vector cImUf = diag(SnCmptMag[faceI] * AImDivInAImU[faceI]) & vector::one;
                        vector cImLf = diag(SnCmptMag[faceI] * AImDivInAImL[faceI]) & vector::one;

                        AImDivInAImD[own[faceI]] += cImUf;
                        AImDivInAImD[nei[faceI]] += cImLf;
                        AImDivInAImU[faceI] -= cImUf;
                        AImDivInAImL[faceI] -= cImLf;
                    }
                }
            }

            // Insert coupling of divA in A
            AEqn.insertBlockCoupling(0, 6, AReDivInARe, true);
            AEqn.insertBlockCoupling(3, 7, AImDivInAIm, true);

            // Assemble coupling of A in divA
            BlockLduSystem<vector, scalar> AReInAReDiv(fvm::UDiv(ARe)); AReInAReDiv *= -rMu0.value();
            BlockLduSystem<vector, scalar> AImInAImDiv(fvm::UDiv(AIm)); AImInAImDiv *= -rMu0.value();

            // Is this necessary?
            // Modify fvm::UDiv() for A*n=0 at conductor interface
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
                CoeffField<vector>::linearTypeField& AReInAReDivD = AReInAReDiv.diag().asLinear();
                CoeffField<vector>::linearTypeField& AReInAReDivU = AReInAReDiv.upper().asLinear();
                CoeffField<vector>::linearTypeField& AReInAReDivL = AReInAReDiv.lower().asLinear();

                CoeffField<vector>::linearTypeField& AImInAImDivD = AImInAImDiv.diag().asLinear();
                CoeffField<vector>::linearTypeField& AImInAImDivU = AImInAImDiv.upper().asLinear();
                CoeffField<vector>::linearTypeField& AImInAImDivL = AImInAImDiv.lower().asLinear();

                // Loop over all faces of conductor interface
                forAll (cnFaces, facei)
                {
                    label faceI = cnFaces[facei];

                    // Eliminate all coeffs
                    // TODO: Optimize
                    {
                        // Real part
                        vector cReUf = diag(SnCmptMag[faceI] * AReInAReDivU[faceI]) & vector::one;
                        vector cReLf = diag(SnCmptMag[faceI] * AReInAReDivL[faceI]) & vector::one;

                        AReInAReDivD[own[faceI]] += cReUf;
                        AReInAReDivD[nei[faceI]] += cReLf;
                        AReInAReDivU[faceI] -= cReUf;
                        AReInAReDivL[faceI] -= cReLf;

                        // Imaginary part
                        vector cImUf = diag(SnCmptMag[faceI] * AImInAImDivU[faceI]) & vector::one;
                        vector cImLf = diag(SnCmptMag[faceI] * AImInAImDivL[faceI]) & vector::one;

                        AImInAImDivD[own[faceI]] += cImUf;
                        AImInAImDivD[nei[faceI]] += cImLf;
                        AImInAImDivU[faceI] -= cImUf;
                        AImInAImDivL[faceI] -= cImLf;
                    }
                }
            }

            // Insert coupling of A in divA
            AEqn.insertBlockCoupling(6, 0, AReInAReDiv, false);
            AEqn.insertBlockCoupling(7, 3, AImInAImDiv, false);
        }
    }



    // Gauge fixing: Solve for gauge fixing potential G
    {
        // Assemble G-equations
        fvScalarMatrix GReEqn
        (
            fvm::laplacian(GRe)
          + fvc::div(ARe)
        );
        fvScalarMatrix GImEqn
        (
            fvm::laplacian(GIm)
          + fvc::div(AIm)
        );

        // Reference settings
        label refCell = mesh.findCell(point(0.00, 0.00, 0.01));
        scalar refValue = 0.0;

        GReEqn.setReference(refCell, refValue);
        GImEqn.setReference(refCell, refValue);

        // Solve G
        GReEqn.solve();
        GImEqn.solve();

        GRe.correctBoundaryConditions();
        GIm.correctBoundaryConditions();

        // TODO: Remove after debug
        volVectorField AReCoulomb("AReCoulomb", ARe + fvc::grad(GRe));
        volVectorField AImCoulomb("AImCoulomb", AIm + fvc::grad(GIm));
        AReCoulomb.write();
        AImCoulomb.write();
    }
