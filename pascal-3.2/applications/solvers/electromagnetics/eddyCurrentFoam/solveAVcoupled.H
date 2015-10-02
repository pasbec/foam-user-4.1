
    {
        // Create scale factor for V-equation compared
        // to A-equation based on dimensional relation
        // dim[A] <=> dim[grad(V)]
        dimensionedScalar scaleV
        (
            "scaleV",
            pow(average(mesh.V()), scalar(1.0/3.0))
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

//         // Coupling terms of electric scalar potential V
//         // in equation for magnetic vector potential A
//         {
//             // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
//
//             // Implicit gradient terms of electric scalar potential V
//             BlockLduSystem<vector, vector> VReInARe(fvm::grad(VRe)); VReInARe *= -sigmaIn;
//             BlockLduSystem<vector, vector> VImInAIm(fvm::grad(VIm)); VImInAIm *= -sigmaIn;
//
// //             // Modify fvm::grad(VRe)/fvm::grad(VIm) to enforce zero normal
// //             // gradient for V (d(V)/dn=0) at conductor interface as inner
// //             // constraint
// //             // TODO: Gradient correction for div(j)=0???
// //             {
// //                 // Mesh addressing
// //                 const unallocLabelList& own = mesh.owner();
// //                 const unallocLabelList& nei = mesh.neighbour();
// //
// //                 // Mesh geometry
// //                 const vectorField& Sf = mesh.Sf();
// //                 const scalarField& magSf = mesh.magSf();
// //
// //                 // Derived mesh geometry
// //                 const vectorField Sn = Sf/magSf;
// //                 const vectorField SnCmptMag = cmptMag(Sn);
// //
// //                 // Reference to parts of blockLduSystem
// //                 CoeffField<vector>::linearTypeField& VReInAReD = VReInARe.diag().asLinear();
// //                 CoeffField<vector>::linearTypeField& VReInAReU = VReInARe.upper().asLinear();
// //                 CoeffField<vector>::linearTypeField& VReInAReL = VReInARe.lower().asLinear();
// //
// //                 CoeffField<vector>::linearTypeField& VImInAImD = VImInAIm.diag().asLinear();
// //                 CoeffField<vector>::linearTypeField& VImInAImU = VImInAIm.upper().asLinear();
// //                 CoeffField<vector>::linearTypeField& VImInAImL = VImInAIm.lower().asLinear();
// //
// //                 // Loop over all faces of conductor interface
// //                 forAll (cnFaces, facei)
// //                 {
// //                     label faceI = cnFaces[facei];
// //
// //                     // Eliminate upper/lower coeffs based on constraint
// //                     // NOTE: This will keep SYMMETRIC matrices and it will
// //                     //       also sligthly improve their diagonal dominance
// //                     {
// //                         // Real part
// //                         vector cReUf = diag(SnCmptMag[faceI] * VReInAReU[faceI]) & vector::one;
// //                         vector cReLf = diag(SnCmptMag[faceI] * VReInAReL[faceI]) & vector::one;
// //
// //                         VReInAReD[own[faceI]] += cReUf;
// //                         VReInAReD[nei[faceI]] += cReLf;
// //                         VReInAReU[faceI] -= cReUf;
// //                         VReInAReL[faceI] -= cReLf;
// //
// //                         // Imaginary part
// //                         vector cImUf = diag(SnCmptMag[faceI] * VImInAImU[faceI]) & vector::one;
// //                         vector cImLf = diag(SnCmptMag[faceI] * VImInAImL[faceI]) & vector::one;
// //
// //                         VImInAImD[own[faceI]] += cImUf;
// //                         VImInAImD[nei[faceI]] += cImLf;
// //                         VImInAImU[faceI] -= cImUf;
// //                         VImInAImL[faceI] -= cImLf;
// //                     }
// //                 }
// //
// //                 // Remove coeffs for non-conducting cells
// //                 // TODO: In how far will this change upper/lower of
// //                 // cells in contact with non-conducting region?
// //                 VReInARe *= indicatorVc;
// //                 VImInAIm *= indicatorVc;
// //             }
//
//             // Insert implicit gradient terms
//             AVEqn.insertBlockCoupling(0, 3, VReInARe, true);
//             AVEqn.insertBlockCoupling(4, 7, VImInAIm, true);
//         }

        // Coupling terms of magnetic vector potential A
        // in equations for electric scalar potential V
        {
            // Implicit gradient constraints
            // at conductor interface
            {
                // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
            }

            // Implicit source terms due to variable sigma
            // with implied condition: div(A)=0
            AVEqn.insertEquationCoupling(3, 4, -scaleV.value() * indicatorVcIn * gradAlphaIn.component(0));
            AVEqn.insertEquationCoupling(3, 5, -scaleV.value() * indicatorVcIn * gradAlphaIn.component(1));
            AVEqn.insertEquationCoupling(3, 6, -scaleV.value() * indicatorVcIn * gradAlphaIn.component(2));
            AVEqn.insertEquationCoupling(7, 0,  scaleV.value() * indicatorVcIn * gradAlphaIn.component(0));
            AVEqn.insertEquationCoupling(7, 1,  scaleV.value() * indicatorVcIn * gradAlphaIn.component(1));
            AVEqn.insertEquationCoupling(7, 2,  scaleV.value() * indicatorVcIn * gradAlphaIn.component(2));
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
