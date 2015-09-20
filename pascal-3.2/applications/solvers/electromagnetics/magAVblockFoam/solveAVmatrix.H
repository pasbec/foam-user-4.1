
    // Create set for conducting cells
    labelHashSet cCellsHashSet;

    forAll(mesh.cells(), cellI)
    {
        if (sigma.internalField()[cellI] > SMALL)
        {
            cCellsHashSet.insert(cellI);
        }
    }

    cellSet cCells
    (
        mesh,
        "cellSet_solver_c",
        cCellsHashSet
    ); cCells.write();

    // Create set for non-conducting cells
    labelHashSet nCellsHashSet;

    forAll(mesh.cells(), cellI)
    {
        if (sigma.internalField()[cellI] < SMALL)
        {
            nCellsHashSet.insert(cellI);
        }
    }

    cellSet nCells
    (
        mesh,
        "cellSet_solver_n",
        nCellsHashSet
    ); nCells.write();

    // Create sets for conducting/non-conducting faces (cn):
    // cnPN) Owner is in conducting region and neighbour
    //       lies in non-conducting region
    // cnNP) Owner is in non-conducting region and
    //       neighbour lies in conducting region
    labelHashSet cnFacesHashSet;
    labelHashSet cnPNfacesHashSet;
    labelHashSet cnNPfacesHashSet;

    {
        // Mesh references
        const unallocLabelList& own = mesh.owner();
        const unallocLabelList& nei = mesh.neighbour();

        // Find and insert faces
        forAll (own, faceI)
        {
            if
            (
                (sigma.internalField()[own[faceI]] > SMALL)
             && (sigma.internalField()[nei[faceI]] < SMALL)
            )
            {
                cnFacesHashSet.insert(faceI);
                cnPNfacesHashSet.insert(faceI);
            }

            if
            (
                (sigma.internalField()[own[faceI]] < SMALL)
             && (sigma.internalField()[nei[faceI]] > SMALL)
            )
            {
                cnFacesHashSet.insert(faceI);
                cnNPfacesHashSet.insert(faceI);
            }
        }
    }

    faceSet cnFaces
    (
        mesh,
        "cellSet_solver_cn",
        cnFacesHashSet
    ); cnFaces.write();

    faceSet cnPNfaces
    (
        mesh,
        "cellSet_solver_cn_PN",
        cnPNfacesHashSet
    ); cnPNfaces.write();

    faceSet cnNPfaces
    (
        mesh,
        "cellSet_solver_cn_NP",
        cnNPfacesHashSet
    ); cnNPfaces.write();



    // Create component fields for A
    volScalarField ARex
    (
        IOobject
        (
            "ARex",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        ),
        mesh,
        dimensioned<scalar>(word(), ARe.dimensions(), pTraits<scalar>::zero),
        fixedValueFvPatchField<scalar>::typeName
    );
    volScalarField ARey
    (
        IOobject
        (
            "ARey",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        ),
        mesh,
        dimensioned<scalar>(word(), ARe.dimensions(), pTraits<scalar>::zero),
        fixedValueFvPatchField<scalar>::typeName
    );
    volScalarField ARez
    (
        IOobject
        (
            "ARez",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        ),
        mesh,
        dimensioned<scalar>(word(), ARe.dimensions(), pTraits<scalar>::zero),
        fixedValueFvPatchField<scalar>::typeName
    );
    volScalarField AImx
    (
        IOobject
        (
            "AImx",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        ),
        mesh,
        dimensioned<scalar>(word(), AIm.dimensions(), pTraits<scalar>::zero),
        fixedValueFvPatchField<scalar>::typeName
    );
    volScalarField AImy
    (
        IOobject
        (
            "AImy",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        ),
        mesh,
        dimensioned<scalar>(word(), AIm.dimensions(), pTraits<scalar>::zero),
        fixedValueFvPatchField<scalar>::typeName
    );
    volScalarField AImz
    (
        IOobject
        (
            "AImz",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        ),
        mesh,
        dimensioned<scalar>(word(), AIm.dimensions(), pTraits<scalar>::zero),
        fixedValueFvPatchField<scalar>::typeName
    );



    // Copy component data
    ARex == ARe.component(0);
    ARey == ARe.component(1);
    ARez == ARe.component(2);
    AImx == AIm.component(0);
    AImy == AIm.component(1);
    AImz == AIm.component(2);

    ARex.correctBoundaryConditions();
    ARey.correctBoundaryConditions();
    ARez.correctBoundaryConditions();
    AImx.correctBoundaryConditions();
    AImy.correctBoundaryConditions();
    AImz.correctBoundaryConditions();



    // Assemble A-eqations
    fvScalarMatrix ARexEqn
    (
        rMu0 * fvm::laplacian(ARex)
//       + alpha * AImx
      + jsRe.component(0)
    );
    fvScalarMatrix AReyEqn
    (
        rMu0 * fvm::laplacian(ARey)
//       + alpha * AImy
      + jsRe.component(1)
    );
    fvScalarMatrix ARezEqn
    (
        rMu0 * fvm::laplacian(ARez)
//       + alpha * AImz
      + jsRe.component(2)
    );
    fvScalarMatrix AImxEqn
    (
        rMu0 * fvm::laplacian(AImx)
//       - alpha * ARex
      + jsIm.component(0)
    );
    fvScalarMatrix AImyEqn
    (
        rMu0 * fvm::laplacian(AImy)
//       - alpha * ARey
      + jsIm.component(1)
    );
    fvScalarMatrix AImzEqn
    (
        rMu0 * fvm::laplacian(AImz)
//       - alpha * ARez
      + jsIm.component(2)
    );



//     // Add explicit terms which are necessary if div(A) is NOT zero
//     volVectorField AReDivGrad("AReDivGrad", fvc::grad(fvc::div(ARe)));
//     volVectorField AImDivGrad("AImDivGrad", fvc::grad(fvc::div(AIm)));
//
//     ARexEqn -= rMu0 * AReDivGrad.component(0);
//     AReyEqn -= rMu0 * AReDivGrad.component(1);
//     ARezEqn -= rMu0 * AReDivGrad.component(2);
//     AImxEqn -= rMu0 * AImDivGrad.component(0);
//     AImyEqn -= rMu0 * AImDivGrad.component(1);
//     AImzEqn -= rMu0 * AImDivGrad.component(2);



//     // Add explicit parts of G
//     volVectorField GReGrad("GReGrad", fvc::grad(GRe));
//     volVectorField GImGrad("GImGrad", fvc::grad(GIm));
//
//     ARexEqn += alpha * GImGrad.component(0);
//     AReyEqn += alpha * GImGrad.component(1);
//     ARezEqn += alpha * GImGrad.component(2);
//     AImxEqn -= alpha * GReGrad.component(0);
//     AImyEqn -= alpha * GReGrad.component(1);
//     AImzEqn -= alpha * GReGrad.component(2);



    {
        // Mesh cell-face addressing
        const unallocLabelList& own = mesh.owner();
        const unallocLabelList& nei = mesh.neighbour();

        // Mesh geometry
        // (Sf/Sn always points from owner to neighbour)
        const vectorField& Cf = mesh.Cf();
        const vectorField& C = mesh.C();
        const vectorField& Sf = mesh.Sf();
        const scalarField& magSf = mesh.magSf();
        const scalarField& d = mesh.deltaCoeffs();

        // Derived mesh geometry
        const vectorField Sn = Sf/magSf;
        const vectorField SnCmptMag = cmptMag(Sn);
        const scalarField SnxMag = mag(Sn.component(0));
        const scalarField SnyMag = mag(Sn.component(1));
        const scalarField SnzMag = mag(Sn.component(2));

        forAll (cnFaces.toc(), facei)
        {
            label faceI = cnFaces.toc()[facei];

            // Calculate Pf and fN
            scalar Pf = mag(Sn[faceI] & (Cf[faceI] - C[own[faceI]]));
            scalar fN = mag(Sn[faceI] & (C[nei[faceI]] - Cf[faceI]));

            // Calculate zero weights
            scalar w = fN / (Pf+fN);
            scalar wOwn = w/(1-w);
            scalar wNei = (1-w)/w;
            scalar wSnx = SnxMag[faceI];
            scalar wSny = SnyMag[faceI];
            scalar wSnz = SnzMag[faceI];

            // Calculate laplacian coeff for this face
            scalar c = rMu0.value() * magSf[faceI] * d[faceI];

            // Version 1: Eliminate upper/neighbour coeffs
            // Note: This will keep a SYMMETRIC matrix and
            //       it will also improve diagonal dominance
            // TODO: Optimize
            ARexEqn.diag()[own[faceI]] -= wSnx * wOwn * c;
            ARexEqn.diag()[nei[faceI]] -= wSnx * wNei * c;
            ARexEqn.upper()[faceI] -= wSnx * c;
            AReyEqn.diag()[own[faceI]] -= wSny * wOwn * c;
            AReyEqn.diag()[nei[faceI]] -= wSny * wNei * c;
            AReyEqn.upper()[faceI] -= wSny * c;
            ARezEqn.diag()[own[faceI]] -= wSnz * wOwn * c;
            ARezEqn.diag()[nei[faceI]] -= wSnz * wNei * c;
            ARezEqn.upper()[faceI] -= wSnz * c;
            AImxEqn.diag()[own[faceI]] -= wSnx * wOwn * c;
            AImxEqn.diag()[nei[faceI]] -= wSnx * wNei * c;
            AImxEqn.upper()[faceI] -= wSnx * c;
            AImyEqn.diag()[own[faceI]] -= wSny * wOwn * c;
            AImyEqn.diag()[nei[faceI]] -= wSny * wNei * c;
            AImyEqn.upper()[faceI] -= wSny * c;
            AImzEqn.diag()[own[faceI]] -= wSnz * wOwn * c;
            AImzEqn.diag()[nei[faceI]] -= wSnz * wNei * c;
            AImzEqn.upper()[faceI] -= wSnz * c;

//             // Version 2: Eliminate diagonal coeffs
//             // Note: This will result in an ASYMMETRIC matrix
//             //       and will dramatically effect the diagonal
//             //       dominance.
//             ARexEqn.diag()[own[faceI]] += wSnx * c;
//             ARexEqn.diag()[nei[faceI]] += wSnx * c;
//             ARexEqn.upper()[faceI] += wSnx * wNei * c;
//             ARexEqn.lower()[faceI] += wSnx * wOwn * c;
//             AReyEqn.diag()[own[faceI]] += wSny * c;
//             AReyEqn.diag()[nei[faceI]] += wSny * c;
//             AReyEqn.upper()[faceI] += wSny * wNei * c;
//             AReyEqn.lower()[faceI] += wSny * wOwn * c;
//             ARezEqn.diag()[own[faceI]] += wSnz * c;
//             ARezEqn.diag()[nei[faceI]] += wSnz * c;
//             ARezEqn.upper()[faceI] += wSnz * wNei * c;
//             ARezEqn.lower()[faceI] += wSnz * wOwn * c;
//             AImxEqn.diag()[own[faceI]] += wSnx * c;
//             AImxEqn.diag()[nei[faceI]] += wSnx * c;
//             AImxEqn.upper()[faceI] += wSnx * wNei * c;
//             AImxEqn.lower()[faceI] += wSnx * wOwn * c;
//             AImyEqn.diag()[own[faceI]] += wSny * c;
//             AImyEqn.diag()[nei[faceI]] += wSny * c;
//             AImyEqn.upper()[faceI] += wSny * wNei * c;
//             AImyEqn.lower()[faceI] += wSny * wOwn * c;
//             AImzEqn.diag()[own[faceI]] += wSnz * c;
//             AImzEqn.diag()[nei[faceI]] += wSnz * c;
//             AImzEqn.upper()[faceI] += wSnz * wNei * c;
//             AImzEqn.lower()[faceI] += wSnz * wOwn * c;

            // TODO: Gradient correction?
            // We have to modifiy th source term, too.
        }
    }

    {

        // References to internal alpha field
        const scalarField& alphaIn = alpha.internalField();

//         // Prepare block system
//         fvBlockMatrix<vector6> AEqn(A);
//
//         // Lookup solver type for AEqn
//         word solverType =
//             mesh.solutionDict().subDict("solvers").
//             subDict(AEqn.psi().name()).lookup("solver");
//
//         // Force diag/upper/lower square for AMG (A is ASYMMETRIC)
//         if (solverType == "AMG")
//         {
//             Field<tensor6>& AEqnD = AEqn.diag().asSquare(); AEqnD = tensor6::zero;
//             Field<tensor6>& AEqnU = AEqn.upper().asSquare(); AEqnU = tensor6::zero;
//             Field<tensor6>& AEqnL = AEqn.lower().asSquare(); AEqnL = tensor6::zero;
//         }

        // Prepare block system for A/divA
        fvBlockMatrix<vector8> AEqn(A);

        {
            // Insert equations for A
            AEqn.insertEquation(0, ARexEqn);
            AEqn.insertEquation(1, AReyEqn);
            AEqn.insertEquation(2, ARezEqn);
            AEqn.insertEquation(3, AImxEqn);
            AEqn.insertEquation(4, AImyEqn);
            AEqn.insertEquation(5, AImzEqn);

            // Add coupling terms for ARe-AIm coupling
            AEqn.insertEquationCoupling(0, 3,  alphaIn);
            AEqn.insertEquationCoupling(3, 0, -alphaIn);
            AEqn.insertEquationCoupling(1, 4,  alphaIn);
            AEqn.insertEquationCoupling(4, 1, -alphaIn);
            AEqn.insertEquationCoupling(2, 5,  alphaIn);
            AEqn.insertEquationCoupling(5, 2, -alphaIn);
        }

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
                fvm::Sp(scalar(1.0),AReDiv)
//               - fvc::div(ARe)
            );

            // Currently only work-around
            // since diagonal insertion
            // in block-matrix does not
            // work correctly
            AReDivEqn.upper() = 0.0;

            fvScalarMatrix AImDivEqn
            (
                fvm::Sp(scalar(1.0),AImDiv)
//               - fvc::div(AIm)
            );

            // Currently only work-around
            // since diagonal insertion
            // in block-matrix does not
            // work correctly
            AImDivEqn.upper() = 0.0;

            // Insert equations for divA
            AEqn.insertEquation(6, AReDivEqn);
            AEqn.insertEquation(7, AImDivEqn);

            // Coupling part of divA in A
            BlockLduSystem<vector, vector> AReDivInARe(fvm::grad(AReDiv)); AReDivInARe *= -rMu0.value();
            BlockLduSystem<vector, vector> AImDivInAIm(fvm::grad(AImDiv)); AImDivInAIm *= -rMu0.value();

            // TODO: Manipulate grad
            {
                // Mesh cell-face addressing
                const unallocLabelList& own = mesh.owner();
                const unallocLabelList& nei = mesh.neighbour();

                // Mesh geometry
                // (Sf/Sn always points from owner to neighbour)
                const vectorField& Sf = mesh.Sf();
                const scalarField& magSf = mesh.magSf();

                // Derived mesh geometry
                const vectorField Sn = Sf/magSf;
                const vectorField SnCmptMag = cmptMag(Sn);
                const scalarField SnxMag = mag(Sn.component(0));
                const scalarField SnyMag = mag(Sn.component(1));
                const scalarField SnzMag = mag(Sn.component(2));

                forAll (cnFaces.toc(), facei)
                {
                    label faceI = cnFaces.toc()[facei];

                    // Calculate zero weights
                    scalar wSnx = SnxMag[faceI];
                    scalar wSny = SnyMag[faceI];
                    scalar wSnz = SnzMag[faceI];

//                     CoeffField<vector>::linearTypeField& d = bs.diag().asLinear();
//                     CoeffField<vector>::linearTypeField& u = bs.upper().asLinear();
//                     CoeffField<vector>::linearTypeField& l = bs.lower().asLinear();

                    // Eliminate all coeffs
                    // TODO: Optimize
                    AReDivInARe.diag().asLinear()[own[faceI]] += AReDivInARe.upper().asLinear()[faceI];
                    AReDivInARe.diag().asLinear()[nei[faceI]] += AReDivInARe.lower().asLinear()[faceI];
                    AReDivInARe.upper().asLinear()[faceI] = vector::zero;
                    AReDivInARe.lower().asLinear()[faceI] = vector::zero;
                    AImDivInAIm.diag().asLinear()[own[faceI]] += AImDivInAIm.upper().asLinear()[faceI];
                    AImDivInAIm.diag().asLinear()[nei[faceI]] += AImDivInAIm.lower().asLinear()[faceI];
                    AImDivInAIm.upper().asLinear()[faceI] = vector::zero;
                    AImDivInAIm.lower().asLinear()[faceI] = vector::zero;

                    // TODO: Direction?
                    vector AReDivInAReUpper = AReDivInARe.upper().asLinear()[faceI];
                    vector AReDivInAReLower = AReDivInARe.lower().asLinear()[faceI];
//                     AReDivInARe.diag().asLinear()[own[faceI]].component(0) += wSnx * vector(1, 0, 0) & AReDivInARe.upper().asLinear()[faceI];
//                     AReDivInARe.diag().asLinear()[own[faceI]].component(1) += wSny * vector(0, 1, 0) & AReDivInARe.upper().asLinear()[faceI];
//                     AReDivInARe.diag().asLinear()[own[faceI]].component(2) += wSnz * vector(0, 0, 1) & AReDivInARe.upper().asLinear()[faceI];
//                     AReDivInARe.diag().asLinear()[nei[faceI]].component(0) += wSnx * vector(1, 0, 0) & AReDivInARe.lower().asLinear()[faceI];
//                     AReDivInARe.diag().asLinear()[nei[faceI]].component(1) += wSny * vector(0, 1, 0) & AReDivInARe.lower().asLinear()[faceI];
//                     AReDivInARe.diag().asLinear()[nei[faceI]].component(2) += wSnz * vector(0, 0, 1) & AReDivInARe.lower().asLinear()[faceI];
//                     AReDivInARe.upper().asLinear()[faceI].component(0) -= wSnx * vector(1, 0, 0) & AReDivInARe.upper().asLinear()[faceI];
//                     AReDivInARe.upper().asLinear()[faceI].component(1) -= wSny * vector(0, 1, 0) & AReDivInARe.upper().asLinear()[faceI];
//                     AReDivInARe.upper().asLinear()[faceI].component(2) -= wSnz * vector(0, 0, 1) & AReDivInARe.upper().asLinear()[faceI];
//                     AReDivInARe.lower().asLinear()[faceI].component(0) -= wSnx * vector(1, 0, 0) & AReDivInARe.lower().asLinear()[faceI];
//                     AReDivInARe.lower().asLinear()[faceI].component(1) -= wSny * vector(0, 1, 0) & AReDivInARe.lower().asLinear()[faceI];
//                     AReDivInARe.lower().asLinear()[faceI].component(2) -= wSnz * vector(0, 0, 1) & AReDivInARe.lower().asLinear()[faceI];
//
//                     AImDivInAIm.diag().asLinear()[own[faceI]].component(0) += wSnx * vector(1, 0, 0) & AImDivInAIm.upper().asLinear()[faceI];
//                     AImDivInAIm.diag().asLinear()[own[faceI]].component(1) += wSny * vector(0, 1, 0) & AImDivInAIm.upper().asLinear()[faceI];
//                     AImDivInAIm.diag().asLinear()[own[faceI]].component(2) += wSnz * vector(0, 0, 1) & AImDivInAIm.upper().asLinear()[faceI];
//                     AImDivInAIm.diag().asLinear()[nei[faceI]].component(0) += wSnx * vector(1, 0, 0) & AImDivInAIm.lower().asLinear()[faceI];
//                     AImDivInAIm.diag().asLinear()[nei[faceI]].component(1) += wSny * vector(0, 1, 0) & AImDivInAIm.lower().asLinear()[faceI];
//                     AImDivInAIm.diag().asLinear()[nei[faceI]].component(2) += wSnz * vector(0, 0, 1) & AImDivInAIm.lower().asLinear()[faceI];
//                     AImDivInAIm.upper().asLinear()[faceI].component(0) -= wSnx * vector(1, 0, 0) & AImDivInAIm.upper().asLinear()[faceI];
//                     AImDivInAIm.upper().asLinear()[faceI].component(1) -= wSny * vector(0, 1, 0) & AImDivInAIm.upper().asLinear()[faceI];
//                     AImDivInAIm.upper().asLinear()[faceI].component(2) -= wSnz * vector(0, 0, 1) & AImDivInAIm.upper().asLinear()[faceI];
//                     AImDivInAIm.lower().asLinear()[faceI].component(0) -= wSnx * vector(1, 0, 0) & AImDivInAIm.lower().asLinear()[faceI];
//                     AImDivInAIm.lower().asLinear()[faceI].component(1) -= wSny * vector(0, 1, 0) & AImDivInAIm.lower().asLinear()[faceI];
//                     AImDivInAIm.lower().asLinear()[faceI].component(2) -= wSnz * vector(0, 0, 1) & AImDivInAIm.lower().asLinear()[faceI];
                }
            }

            AEqn.insertBlockCoupling(0, 6, AReDivInARe, true);
            AEqn.insertBlockCoupling(3, 7, AImDivInAIm, true);

            // Coupling part of A in divA
            BlockLduSystem<vector, scalar> AReInAReDiv(fvm::UDiv(ARe)); AReInAReDiv *= -1.0;
            BlockLduSystem<vector, scalar> AImInAImDiv(fvm::UDiv(AIm)); AImInAImDiv *= -1.0;

            // TODO: Manipulate UDiv
            {
            }

            AEqn.insertBlockCoupling(6, 0, AReInAReDiv, false);
            AEqn.insertBlockCoupling(7, 3, AImInAImDiv, false);
        }

        // Block coupled solver call to solve A
        AEqn.solve();

        // Retrieve solution
        AEqn.retrieveSolution(0, ARex.internalField());
        AEqn.retrieveSolution(1, ARey.internalField());
        AEqn.retrieveSolution(2, ARez.internalField());
        AEqn.retrieveSolution(3, AImx.internalField());
        AEqn.retrieveSolution(4, AImy.internalField());
        AEqn.retrieveSolution(5, AImz.internalField());

        // TODO: Remove after debug
        AEqn.retrieveSolution(6, AReDiv.internalField());
        AEqn.retrieveSolution(7, AImDiv.internalField());
        AReDiv.write();
        AImDiv.write();

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

        // TODO: This has to be the same as the
        // component boundary conditions otherwise
        // it gives non-sense. This means we currently
        // cannot use directionMixed bc!!!
        ARe.correctBoundaryConditions();
        AIm.correctBoundaryConditions();
    }

//     {
//         // Assemble G-eqations
//         fvScalarMatrix GReEqn(fvm::laplacian(GRe) + fvc::div(ARe));
//         fvScalarMatrix GImEqn(fvm::laplacian(GIm) + fvc::div(AIm));
//
//         // Solve G
//         GReEqn.solve();
//         GImEqn.solve();
//
//         GRe.correctBoundaryConditions();
//         GIm.correctBoundaryConditions();
//
//         // Correct A
//         AIm += alpha * fvc::grad(GIm);
//         ARe -= alpha * fvc::grad(GRe);
//     }

    jRe ==   alpha * AIm;
    jIm == - alpha * ARe;

    BRe == fvc::curl(ARe);
    BIm == fvc::curl(AIm);
