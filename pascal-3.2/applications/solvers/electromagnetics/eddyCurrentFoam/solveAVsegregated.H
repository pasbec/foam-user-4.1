
    // Solve A
    {
        // Coupling terms of electric scalar potential V
        // in equation for magnetic vector potential A
        {
            // Explicit gradient terms of electric scalar potential V
            AReEqn -= sigma * VReGrad;
            AImEqn -= sigma * VImGrad;
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
            AVEqn.insertEquationCoupling(0, 4,  alphaIn);
            AVEqn.insertEquationCoupling(4, 0, -alphaIn);
            AVEqn.insertEquationCoupling(1, 5,  alphaIn);
            AVEqn.insertEquationCoupling(5, 1, -alphaIn);
            AVEqn.insertEquationCoupling(2, 6,  alphaIn);
            AVEqn.insertEquationCoupling(6, 2, -alphaIn);
        }

        // Solve A
        AEqn.solve();

        // Retrieve solution
        AEqn.retrieveSolution(0, ARe.internalField());
        AEqn.retrieveSolution(3, AIm.internalField());

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
                const vectorField& Cf = mesh.Cf();
                const vectorField& C = mesh.C();
                const vectorField& Sf = mesh.Sf();
                const scalarField& magSf = mesh.magSf();
                const vectorField Sn = Sf/magSf;

                // ...
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
                    scalar cP = omega.value() * sSigmaV[faceI] * magSf[faceI] * wP;
                    scalar cN = omega.value() * sSigmaV[faceI] * magSf[faceI] * wN;

                    VReEqn.source()[own[faceI]] -= cP * (Sn[faceI] & AIm[own[faceI]]);
                    VReEqn.source()[own[faceI]] -= cN * (Sn[faceI] & AIm[nei[faceI]]);
                    VReEqn.source()[nei[faceI]] += cP * (Sn[faceI] & AIm[own[faceI]]);
                    VReEqn.source()[nei[faceI]] += cN * (Sn[faceI] & AIm[nei[faceI]]);

                    VImEqn.source()[own[faceI]] += cP * (Sn[faceI] & ARe[own[faceI]]);
                    VImEqn.source()[own[faceI]] += cN * (Sn[faceI] & ARe[nei[faceI]]);
                    VImEqn.source()[nei[faceI]] -= cP * (Sn[faceI] & ARe[own[faceI]]);
                    VImEqn.source()[nei[faceI]] -= cN * (Sn[faceI] & ARe[nei[faceI]]);
                }
            }

            // Explicit source terms due to variable sigma
            // (force div(A)=0)
            VReEqn -= indicatorVc * omega * (gradAlpha & AIm);
            VImEqn += indicatorVc * omega * (gradAlpha & ARe);
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
