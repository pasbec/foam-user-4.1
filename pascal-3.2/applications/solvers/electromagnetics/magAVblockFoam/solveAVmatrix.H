
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



    // Assemble A-equations
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



    // TODO: Remove after debug
    volScalarField AReDiv("AReDiv", fvc::div(ARe));
    volScalarField AImDiv("AImDiv", fvc::div(AIm));
    AReDiv.write();
    AImDiv.write();

    // Add explicit terms which are necessary if div(A) is NOT zero
    volVectorField AReDivGrad("AReDivGrad", fvc::grad(fvc::div(ARe)));
    volVectorField AImDivGrad("AImDivGrad", fvc::grad(fvc::div(AIm)));

    ARexEqn -= rMu0 * AReDivGrad.component(0);
    AReyEqn -= rMu0 * AReDivGrad.component(1);
    ARezEqn -= rMu0 * AReDivGrad.component(2);
    AImxEqn -= rMu0 * AImDivGrad.component(0);
    AImyEqn -= rMu0 * AImDivGrad.component(1);
    AImzEqn -= rMu0 * AImDivGrad.component(2);

    // TODO: Remove after debug
    AReDivGrad.write();
    AImDivGrad.write();



//     // Add explicit parts of electric scalar potential V
//     ARexEqn -= sigma * VReGrad.component(0);
//     AReyEqn -= sigma * VReGrad.component(1);
//     ARezEqn -= sigma * VReGrad.component(2);
//     AImxEqn -= sigma * VImGrad.component(0);
//     AImyEqn -= sigma * VImGrad.component(1);
//     AImzEqn -= sigma * VImGrad.component(2);



//     // Add explicit parts of gauge fixing potential G
//     ARexEqn += alpha * GImGrad.component(0);
//     AReyEqn += alpha * GImGrad.component(1);
//     ARezEqn += alpha * GImGrad.component(2);
//     AImxEqn -= alpha * GReGrad.component(0);
//     AImyEqn -= alpha * GReGrad.component(1);
//     AImzEqn -= alpha * GReGrad.component(2);
//
//     ARexEqn += rMu0 * fvc::laplacian(GReGrad.component(0), "laplacian(ARex)");
//     AReyEqn += rMu0 * fvc::laplacian(GReGrad.component(1), "laplacian(ARey)");
//     ARezEqn += rMu0 * fvc::laplacian(GReGrad.component(2), "laplacian(ARez)");
//     AImxEqn += rMu0 * fvc::laplacian(GImGrad.component(0), "laplacian(AImx)");
//     AImyEqn += rMu0 * fvc::laplacian(GImGrad.component(1), "laplacian(AImy)");
//     AImzEqn += rMu0 * fvc::laplacian(GImGrad.component(2), "laplacian(AImz)");
//
// // //     ARexEqn -= rMu0 * fvc::grad(fvc::laplacian(GRe))().component(0);
// // //     AReyEqn -= rMu0 * fvc::grad(fvc::laplacian(GRe))().component(1);
// // //     ARezEqn -= rMu0 * fvc::grad(fvc::laplacian(GRe))().component(2);
// // //     AImxEqn -= rMu0 * fvc::grad(fvc::laplacian(GIm))().component(0);
// // //     AImyEqn -= rMu0 * fvc::grad(fvc::laplacian(GIm))().component(1);
// // //     AImzEqn -= rMu0 * fvc::grad(fvc::laplacian(GIm))().component(2);
// //
// //     ARexEqn -= rMu0 * fvc::grad(fvc::div(ARe + GReGrad, "div(ARe)"))().component(0);
// //     AReyEqn -= rMu0 * fvc::grad(fvc::div(ARe + GReGrad, "div(ARe)"))().component(1);
// //     ARezEqn -= rMu0 * fvc::grad(fvc::div(ARe + GReGrad, "div(ARe)"))().component(2);
// //     AImxEqn -= rMu0 * fvc::grad(fvc::div(AIm + GImGrad, "div(AIm)"))().component(0);
// //     AImyEqn -= rMu0 * fvc::grad(fvc::div(AIm + GImGrad, "div(AIm)"))().component(1);
// //     AImzEqn -= rMu0 * fvc::grad(fvc::div(AIm + GImGrad, "div(AIm)"))().component(2);
//
//     ARe -= GReGrad;
//     AIm -= GImGrad;



    // Modify fvm::laplacian() for A*n=0 at c/n-interface
    {
        // Mesh cell-face addressing
        const unallocLabelList& own = mesh.owner();
        const unallocLabelList& nei = mesh.neighbour();

        // Mesh geometry
        const vectorField& Cf = mesh.Cf();
        const vectorField& C = mesh.C();
        const vectorField& Sf = mesh.Sf();
        const scalarField& magSf = mesh.magSf();

        // Derived mesh geometry
        const vectorField Sn = Sf/magSf;
        const vectorField SnCmptMag = cmptMag(Sn);
        const scalarField SnxMag = mag(Sn.component(0));
        const scalarField SnyMag = mag(Sn.component(1));
        const scalarField SnzMag = mag(Sn.component(2));

        // Loop over all faces of c/n-interface
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

            // Eliminate upper/neighbour coeffs
            // Note: This will keep a SYMMETRIC matrix and
            //       it will also improve diagonal dominance
            // TODO: Optimize
            scalar wSnARexEqnUf = wSnx * ARexEqn.upper()[faceI];
            ARexEqn.diag()[own[faceI]] -= wSnARexEqnUf * wOwn;
            ARexEqn.diag()[nei[faceI]] -= wSnARexEqnUf * wNei;
            ARexEqn.upper()[faceI] -= wSnARexEqnUf;
            scalar wSnAReyEqnUf = wSny * AReyEqn.upper()[faceI];
            AReyEqn.diag()[own[faceI]] -= wSnAReyEqnUf * wOwn;
            AReyEqn.diag()[nei[faceI]] -= wSnAReyEqnUf * wNei;
            AReyEqn.upper()[faceI] -= wSnAReyEqnUf;
            scalar wSnARezEqnUf = wSnz * AReyEqn.upper()[faceI];
            ARezEqn.diag()[own[faceI]] -= wSnARezEqnUf * wOwn;
            ARezEqn.diag()[nei[faceI]] -= wSnARezEqnUf * wNei;
            ARezEqn.upper()[faceI] -= wSnARezEqnUf;
            scalar wSnAImxEqnUf = wSnx * AImxEqn.upper()[faceI];
            AImxEqn.diag()[own[faceI]] -= wSnAImxEqnUf * wOwn;
            AImxEqn.diag()[nei[faceI]] -= wSnAImxEqnUf * wNei;
            AImxEqn.upper()[faceI] -= wSnAImxEqnUf;
            scalar wSnAImyEqnUf = wSny * AImyEqn.upper()[faceI];
            AImyEqn.diag()[own[faceI]] -= wSnAImyEqnUf * wOwn;
            AImyEqn.diag()[nei[faceI]] -= wSnAImyEqnUf * wNei;
            AImyEqn.upper()[faceI] -= wSnAImyEqnUf;
            scalar wSnAImzEqnUf = wSnz * AImyEqn.upper()[faceI];
            AImzEqn.diag()[own[faceI]] -= wSnAImzEqnUf * wOwn;
            AImzEqn.diag()[nei[faceI]] -= wSnAImzEqnUf * wNei;
            AImzEqn.upper()[faceI] -= wSnAImzEqnUf;

            // TODO: What about gradient correction?
            // We have to modifiy the source term, too.
        }
    }

    // Solve A
    {

        // References to internal alpha field
        const scalarField& alphaIn = alpha.internalField();

        // Prepare block system for A/divA
        fvBlockMatrix<vector6> AEqn(A);
//         fvBlockMatrix<vector8> AEqn(A);

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

//         volScalarField AReDiv
//         (
//             IOobject
//             (
//                 "AReDiv",
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE
//             ),
//             fvc::div(ARe),
//             fixedValueFvPatchField<scalar>::typeName
//         );
//
//         volScalarField AImDiv
//         (
//             IOobject
//             (
//                 "AImDiv",
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE
//             ),
//             fvc::div(AIm),
//             fixedValueFvPatchField<scalar>::typeName
//         );
//
//         {
//             fvScalarMatrix AReDivEqn
//             (
//                 rMu0 * fvm::Sp(scalar(1.0),AReDiv)
// //            - rMu0 * fvc::div(ARe)
//             );
//
//             // Currently only work-around
//             // since diagonal insertion
//             // in block-matrix does not
//             // work correctly
//             AReDivEqn.upper() = 0.0;
//
//             fvScalarMatrix AImDivEqn
//             (
//                 rMu0 * fvm::Sp(scalar(1.0),AImDiv)
// //            - rMu0 * fvc::div(AIm)
//             );
//
//             // Currently only work-around
//             // since diagonal insertion
//             // in block-matrix does not
//             // work correctly
//             AImDivEqn.upper() = 0.0;
//
//             // Insert equations for divA
//             AEqn.insertEquation(6, AReDivEqn);
//             AEqn.insertEquation(7, AImDivEqn);
//
//             // Assemble coupling of divA in A
//             BlockLduSystem<vector, vector> AReDivInARe(fvm::grad(AReDiv)); AReDivInARe *= -rMu0.value();
//             BlockLduSystem<vector, vector> AImDivInAIm(fvm::grad(AImDiv)); AImDivInAIm *= -rMu0.value();
//
//             // Modify fvm::grad() for A*n=0 at c/n-interface
//             {
//                 // Mesh addressing
//                 const unallocLabelList& own = mesh.owner();
//                 const unallocLabelList& nei = mesh.neighbour();
//
//                 // Mesh geometry
//                 const vectorField& Sf = mesh.Sf();
//                 const scalarField& magSf = mesh.magSf();
//
//                 // Derived mesh geometry
//                 const vectorField Sn = Sf/magSf;
//                 const vectorField SnCmptMag = cmptMag(Sn);
//
//                 // Reference to parts of blockLduSystem
//                 CoeffField<vector>::linearTypeField& AReDivInAReD = AReDivInARe.diag().asLinear();
//                 CoeffField<vector>::linearTypeField& AReDivInAReU = AReDivInARe.upper().asLinear();
//                 CoeffField<vector>::linearTypeField& AReDivInAReL = AReDivInARe.lower().asLinear();
//
//                 CoeffField<vector>::linearTypeField& AImDivInAImD = AImDivInAIm.diag().asLinear();
//                 CoeffField<vector>::linearTypeField& AImDivInAImU = AImDivInAIm.upper().asLinear();
//                 CoeffField<vector>::linearTypeField& AImDivInAImL = AImDivInAIm.lower().asLinear();
//
//                 // Loop over all faces of c/n-interface
//                 forAll (cnFaces.toc(), facei)
//                 {
//                     label faceI = cnFaces.toc()[facei];
//
//                     // Eliminate all coeffs
//                     // TODO: Optimize
//                     {
//                         // Real part
//                         vector cReUf = diag(SnCmptMag[faceI] * AReDivInAReU[faceI]) & vector::one;
//                         vector cReLf = diag(SnCmptMag[faceI] * AReDivInAReL[faceI]) & vector::one;
//
//                         AReDivInAReD[own[faceI]] += cReUf;
//                         AReDivInAReD[nei[faceI]] += cReLf;
//                         AReDivInAReU[faceI] -= cReUf;
//                         AReDivInAReL[faceI] -= cReLf;
//
//                         // Imaginary part
//                         vector cImUf = diag(SnCmptMag[faceI] * AImDivInAImU[faceI]) & vector::one;
//                         vector cImLf = diag(SnCmptMag[faceI] * AImDivInAImL[faceI]) & vector::one;
//
//                         AImDivInAImD[own[faceI]] += cImUf;
//                         AImDivInAImD[nei[faceI]] += cImLf;
//                         AImDivInAImU[faceI] -= cImUf;
//                         AImDivInAImL[faceI] -= cImLf;
//                     }
//                 }
//             }
//
//             // Insert coupling of divA in A
//             AEqn.insertBlockCoupling(0, 6, AReDivInARe, true);
//             AEqn.insertBlockCoupling(3, 7, AImDivInAIm, true);
//
//             // Assemble coupling of A in divA
//             BlockLduSystem<vector, scalar> AReInAReDiv(fvm::UDiv(ARe)); AReInAReDiv *= -rMu0.value();
//             BlockLduSystem<vector, scalar> AImInAImDiv(fvm::UDiv(AIm)); AImInAImDiv *= -rMu0.value();
//
//             // Is this necessary?
//             // Modify fvm::UDiv() for A*n=0 at c/n-interface
//             {
//                 // Mesh addressing
//                 const unallocLabelList& own = mesh.owner();
//                 const unallocLabelList& nei = mesh.neighbour();
//
//                 // Mesh geometry
//                 const vectorField& Sf = mesh.Sf();
//                 const scalarField& magSf = mesh.magSf();
//
//                 // Derived mesh geometry
//                 const vectorField Sn = Sf/magSf;
//                 const vectorField SnCmptMag = cmptMag(Sn);
//
//                 // Reference to parts of blockLduSystem
//                 CoeffField<vector>::linearTypeField& AReInAReDivD = AReInAReDiv.diag().asLinear();
//                 CoeffField<vector>::linearTypeField& AReInAReDivU = AReInAReDiv.upper().asLinear();
//                 CoeffField<vector>::linearTypeField& AReInAReDivL = AReInAReDiv.lower().asLinear();
//
//                 CoeffField<vector>::linearTypeField& AImInAImDivD = AImInAImDiv.diag().asLinear();
//                 CoeffField<vector>::linearTypeField& AImInAImDivU = AImInAImDiv.upper().asLinear();
//                 CoeffField<vector>::linearTypeField& AImInAImDivL = AImInAImDiv.lower().asLinear();
//
//                 // Loop over all faces of c/n-interface
//                 forAll (cnFaces.toc(), facei)
//                 {
//                     label faceI = cnFaces.toc()[facei];
//
//                     // Eliminate all coeffs
//                     // TODO: Optimize
//                     {
//                         // Real part
//                         vector cReUf = diag(SnCmptMag[faceI] * AReInAReDivU[faceI]) & vector::one;
//                         vector cReLf = diag(SnCmptMag[faceI] * AReInAReDivL[faceI]) & vector::one;
//
//                         AReInAReDivD[own[faceI]] += cReUf;
//                         AReInAReDivD[nei[faceI]] += cReLf;
//                         AReInAReDivU[faceI] -= cReUf;
//                         AReInAReDivL[faceI] -= cReLf;
//
//                         // Imaginary part
//                         vector cImUf = diag(SnCmptMag[faceI] * AImInAImDivU[faceI]) & vector::one;
//                         vector cImLf = diag(SnCmptMag[faceI] * AImInAImDivL[faceI]) & vector::one;
//
//                         AImInAImDivD[own[faceI]] += cImUf;
//                         AImInAImDivD[nei[faceI]] += cImLf;
//                         AImInAImDivU[faceI] -= cImUf;
//                         AImInAImDivL[faceI] -= cImLf;
//                     }
//                 }
//             }
//
//             // Insert coupling of A in divA
//             AEqn.insertBlockCoupling(6, 0, AReInAReDiv, false);
//             AEqn.insertBlockCoupling(7, 3, AImInAImDiv, false);
//         }

        // Block coupled solver call to solve A
        Info << "BlockMatrix solver call for A" << endl;
        AEqn.solve();

        // Retrieve solution
        AEqn.retrieveSolution(0, ARex.internalField());
        AEqn.retrieveSolution(1, ARey.internalField());
        AEqn.retrieveSolution(2, ARez.internalField());
        AEqn.retrieveSolution(3, AImx.internalField());
        AEqn.retrieveSolution(4, AImy.internalField());
        AEqn.retrieveSolution(5, AImz.internalField());

//         // TODO: Remove after debug
//         AEqn.retrieveSolution(6, AReDiv.internalField());
//         AEqn.retrieveSolution(7, AImDiv.internalField());
//         AReDiv.write();
//         AImDiv.write();

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

//     // Solve G in non-conducting subregion
//     {
//         // Create submesh
//         fvMeshSubset cMeshSubSet
//         (
//             IOobject
//             (
//                 "conductor",
//                 mesh.time().timeName(),
//                 mesh,
//                 IOobject::NO_READ
//             ),
//             mesh
//         );
//
//         // Initialize large cell subset
//         cMeshSubSet.setLargeCellSubset(cCellsHashSet);
//
//         // Create sub-mesh
//         fvMesh& cMesh = cMeshSubSet.subMesh();
//
//         // References for maps
//         const labelList& cCellMap = cMeshSubSet.cellMap();
//
//         volVectorField cARe
//         (
//             IOobject
//             (
//                 "ARe",
//                 cMesh.time().timeName(),
//                 cMesh,
//                 IOobject::NO_READ
//             ),
//             cMesh,
//             dimensioned<vector>(word(), ARe.dimensions(), pTraits<vector>::zero),
//             zeroGradientFvPatchField<vector>::typeName
//         );
//
//         volVectorField cAIm
//         (
//             IOobject
//             (
//                 "AIm",
//                 cMesh.time().timeName(),
//                 cMesh,
//                 IOobject::NO_READ
//             ),
//             cMesh,
//             dimensionedVector(word(), AIm.dimensions(), vector::zero),
//             zeroGradientFvPatchField<vector>::typeName
//         );
//
//         volScalarField cGRe
//         (
//             IOobject
//             (
//                 "GRe",
//                 cMesh.time().timeName(),
//                 cMesh,
//                 IOobject::NO_READ
//             ),
//             cMesh,
//             dimensionedScalar(word(), GRe.dimensions(), pTraits<scalar>::zero),
//             zeroGradientFvPatchField<scalar>::typeName
//         );
//
//         volScalarField cGIm
//         (
//             IOobject
//             (
//                 "GIm",
//                 cMesh.time().timeName(),
//                 cMesh,
//                 IOobject::NO_READ
//             ),
//             cMesh,
//             dimensionedScalar(word(), GIm.dimensions(), pTraits<scalar>::zero),
//             zeroGradientFvPatchField<scalar>::typeName
//         );
//
//         // Map internal values of volume fields
//         // and correct their boundary values
//         {
//             forAll (cMesh.C(), cellI)
//             {
//                 cARe.internalField()[cellI] = ARe.internalField()[cCellMap[cellI]];
//                 cAIm.internalField()[cellI] = AIm.internalField()[cCellMap[cellI]];
//
//                 cGRe.internalField()[cellI] = GRe.internalField()[cCellMap[cellI]];
//                 cGIm.internalField()[cellI] = GIm.internalField()[cCellMap[cellI]];
//             }
//
//             cARe.correctBoundaryConditions();
//             cAIm.correctBoundaryConditions();
//
//             cGRe.correctBoundaryConditions();
//             cGIm.correctBoundaryConditions();
//         }
//
//         // Solve G equation
//         {
//             // Assemble V-equations
//             fvScalarMatrix cGReEqn(fvm::laplacian(cGRe) + fvc::div(cARe));
//             fvScalarMatrix cGImEqn(fvm::laplacian(cGIm) + fvc::div(cAIm));
//
//             label refCell = cMesh.findCell(point(0.00, 0.00, 0.01));
//             scalar refValue = 0.0;
//
//             cGReEqn.setReference(refCell, refValue);
//             cGImEqn.setReference(refCell, refValue);
//
//             // Solve V
//             cGReEqn.solve();
//             cGImEqn.solve();
//
//             cGRe.correctBoundaryConditions();
//             cGIm.correctBoundaryConditions();
//         }
//
//         // Reverse map fields
//         {
//             // Calculate gradient of electric scalar potential
//             volVectorField cGReGrad("GReGrad", fvc::grad(cGRe));
//             volVectorField cGImGrad("GImGrad", fvc::grad(cGIm));
//
//             // Correct A with gradient of gauge potential G
//             cARe += cGReGrad;
//             cAIm += cGImGrad;
//
//             forAll (cMesh.C(), cellI)
//             {
//                 GRe[cCellMap[cellI]] = cGRe[cellI];
//                 GIm[cCellMap[cellI]] = cGIm[cellI];
//
//                 GReGrad[cCellMap[cellI]] = cGReGrad[cellI];
//                 GImGrad[cCellMap[cellI]] = cGImGrad[cellI];
//
//                 ARe[cCellMap[cellI]] = cARe[cellI];
//                 AIm[cCellMap[cellI]] = cAIm[cellI];
//             }
//         }
//     }
//
//     // Solve G in non-conducting subregion
//     {
//         // Create submesh
//         fvMeshSubset nMeshSubSet
//         (
//             IOobject
//             (
//                 "air",
//                 mesh.time().timeName(),
//                 mesh,
//                 IOobject::NO_READ
//             ),
//             mesh
//         );
//
//         // Initialize large cell subset
//         nMeshSubSet.setLargeCellSubset(nCellsHashSet);
//
//         // Create sub-mesh
//         fvMesh& nMesh = nMeshSubSet.subMesh();
//
//         // References for maps
//         const labelList& nCellMap = nMeshSubSet.cellMap();
//
//         volVectorField nARe
//         (
//             IOobject
//             (
//                 "ARe",
//                 nMesh.time().timeName(),
//                 nMesh,
//                 IOobject::NO_READ
//             ),
//             nMesh,
//             dimensioned<vector>(word(), ARe.dimensions(), pTraits<vector>::zero),
//             zeroGradientFvPatchField<vector>::typeName
//         );
//
//         volVectorField nAIm
//         (
//             IOobject
//             (
//                 "AIm",
//                 nMesh.time().timeName(),
//                 nMesh,
//                 IOobject::NO_READ
//             ),
//             nMesh,
//             dimensionedVector(word(), AIm.dimensions(), vector::zero),
//             zeroGradientFvPatchField<vector>::typeName
//         );
//
//         volScalarField nGRe
//         (
//             IOobject
//             (
//                 "GRe",
//                 nMesh.time().timeName(),
//                 nMesh,
//                 IOobject::NO_READ
//             ),
//             nMesh,
//             dimensionedScalar(word(), GRe.dimensions(), pTraits<scalar>::zero),
//             zeroGradientFvPatchField<scalar>::typeName
//         );
//
//         volScalarField nGIm
//         (
//             IOobject
//             (
//                 "GIm",
//                 nMesh.time().timeName(),
//                 nMesh,
//                 IOobject::NO_READ
//             ),
//             nMesh,
//             dimensionedScalar(word(), GIm.dimensions(), pTraits<scalar>::zero),
//             zeroGradientFvPatchField<scalar>::typeName
//         );
//
//         // Map internal values of volume fields
//         // and correct their boundary values
//         {
//             forAll (nMesh.C(), cellI)
//             {
//                 nARe.internalField()[cellI] = ARe.internalField()[nCellMap[cellI]];
//                 nAIm.internalField()[cellI] = AIm.internalField()[nCellMap[cellI]];
//
//                 nGRe.internalField()[cellI] = GRe.internalField()[nCellMap[cellI]];
//                 nGIm.internalField()[cellI] = GIm.internalField()[nCellMap[cellI]];
//             }
//
//             nARe.correctBoundaryConditions();
//             nAIm.correctBoundaryConditions();
//
//             nGRe.correctBoundaryConditions();
//             nGIm.correctBoundaryConditions();
//         }
//
//         // Solve G equation
//         {
//             // Assemble V-equations
//             fvScalarMatrix nGReEqn(fvm::laplacian(nGRe) + fvc::div(nARe));
//             fvScalarMatrix nGImEqn(fvm::laplacian(nGIm) + fvc::div(nAIm));
//
//             label refCell = nMesh.findCell(point(0.15, 0.00, 0.01));
//             scalar refValue = 0.0;
//
//             nGReEqn.setReference(refCell, refValue);
//             nGImEqn.setReference(refCell, refValue);
//
//             // Solve V
//             nGReEqn.solve();
//             nGImEqn.solve();
//
//             nGRe.correctBoundaryConditions();
//             nGIm.correctBoundaryConditions();
//         }
//
//         // Reverse map fields
//         {
//             // Calculate gradient of electric scalar potential
//             volVectorField nGReGrad("GReGrad", fvc::grad(nGRe));
//             volVectorField nGImGrad("GImGrad", fvc::grad(nGIm));
//
//             // Correct A with gradient of gauge potential G
//             nARe += nGReGrad;
//             nAIm += nGImGrad;
//
//             forAll (nMesh.C(), cellI)
//             {
//                 GRe[nCellMap[cellI]] = nGRe[cellI];
//                 GIm[nCellMap[cellI]] = nGIm[cellI];
//
//                 GReGrad[nCellMap[cellI]] = nGReGrad[cellI];
//                 GImGrad[nCellMap[cellI]] = nGImGrad[cellI];
//
//                 ARe[nCellMap[cellI]] = nARe[cellI];
//                 AIm[nCellMap[cellI]] = nAIm[cellI];
//             }
//         }
//     }

    // Gauge fixing: Solve for gauge fixing potential G
    {
//         surfaceScalarField sARe("sARe", fvc::interpolate(ARe) & mesh.Sf());
//         surfaceScalarField sAIm("sAIm", fvc::interpolate(AIm) & mesh.Sf());
//
//         // Modify surface A for A*n=0 at c/n-interface
//         {
//             // Reference to internal field
//             scalarField& sAReIn = sARe.internalField();
//             scalarField& sAImIn = sAIm.internalField();
//
//             // Loop over all faces of c/n-interface
//             forAll (cnFaces.toc(), facei)
//             {
//                 label faceI = cnFaces.toc()[facei];
//
//                 {
//                     // Set face values to zero
//                     sAReIn[faceI] = 0.0;
//                     sAImIn[faceI] = 0.0;
//                 }
//             }
//         }
//
//         // Calculate corrected div(A)
//         volScalarField AReDivG("AReDivG", fvc::div(sARe));
//         volScalarField AImDivG("AImDivG", fvc::div(sAIm));
//
// //         // Modify div(A) for grad(G)=0 at c/n-interface
// //         {
// //             // Mesh cell-face addressing
// //             const unallocLabelList& own = mesh.owner();
// //             const unallocLabelList& nei = mesh.neighbour();
// //
// //             // Mesh geometry
// //             const vectorField& Sf = mesh.Sf();
// //
// //             // Reference to internal field
// //             scalarField& AReDivGin = AReDivG.internalField();
// //             scalarField& AImDivGin = AImDivG.internalField();
// //
// //             // Loop over all faces of c/n-interface
// //             forAll (cnFaces.toc(), facei)
// //             {
// //                 label faceI = cnFaces.toc()[facei];
// //
// //                 {
// //                     // Remove summands based on surface values
// //                     AReDivGin[own[faceI]] -= sARe[faceI];
// //                     AReDivGin[nei[faceI]] += sARe[faceI];
// //                     AImDivGin[own[faceI]] -= sAIm[faceI];
// //                     AImDivGin[nei[faceI]] += sAIm[faceI];
// //
// //                     // Add summands based on cell values
// //                     AReDivGin[own[faceI]] += ARe[own[faceI]] & Sf[faceI];
// //                     AReDivGin[nei[faceI]] -= ARe[nei[faceI]] & Sf[faceI];
// //                     AImDivGin[own[faceI]] += AIm[own[faceI]] & Sf[faceI];
// //                     AImDivGin[nei[faceI]] -= AIm[nei[faceI]] & Sf[faceI];
// //                 }
// //             }
// //         }
//
//         // TODO: Remove after debug
//         AReDivG.write();
//         AImDivG.write();
//
//         // Assemble G-equations
//         fvScalarMatrix GReEqn(fvm::laplacian(GRe) + AReDivG);
//         fvScalarMatrix GImEqn(fvm::laplacian(GIm) + AImDivG);

        // Assemble G-equations
        fvScalarMatrix GReEqn(fvm::laplacian(GRe) + fvc::div(ARe));
        fvScalarMatrix GImEqn(fvm::laplacian(GIm) + fvc::div(AIm));

//         // Modify fvm::laplacian() to get grad(G)=0 at c/n-interface
//         {
//             // Mesh cell-face addressing
//             const unallocLabelList& own = mesh.owner();
//             const unallocLabelList& nei = mesh.neighbour();
//
//             // Loop over all faces of c/n-interface
//             forAll (cnFaces.toc(), facei)
//             {
//                 label faceI = cnFaces.toc()[facei];
//
//                 // Eliminate upper/neighbour coeffs
//                 // Note: This will keep a SYMMETRIC matrix and
//                 //       it will also keep diagonal dominance
//                 scalar GReEqnUf = GReEqn.upper()[faceI];
//                 GReEqn.diag()[own[faceI]] += GReEqnUf;
//                 GReEqn.diag()[nei[faceI]] += GReEqnUf;
//                 GReEqn.upper()[faceI] -= GReEqnUf;
//                 scalar GImEqnUf = GImEqn.upper()[faceI];
//                 GImEqn.diag()[own[faceI]] += GImEqnUf;
//                 GImEqn.diag()[nei[faceI]] += GImEqnUf;
//                 GImEqn.upper()[faceI] -= GImEqnUf;
//             }
//         }

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

        GReGrad = fvc::grad(GRe);
        GImGrad = fvc::grad(GIm);

        // Correct A with gradient of gauge potential G
//         ARe += GReGrad;
//         AIm += GImGrad;
        
        // TODO: Remove after debug
        volVectorField AReCoulomb("AReCoulomb", ARe + GReGrad);
        volVectorField AImCoulomb("AImCoulomb", AIm + GImGrad);
        AReCoulomb.write();
        AImCoulomb.write();
    }

//     // Solve V in conducting subregion
//     {
//         // Create submesh
//         fvMeshSubset cMeshSubSet
//         (
//             IOobject
//             (
//                 "conductor",
//                 mesh.time().timeName(),
//                 mesh,
//                 IOobject::NO_READ
//             ),
//             mesh
//         );
//
//         // Initialize large cell subset
//         cMeshSubSet.setLargeCellSubset(cCellsHashSet);
//
//         // Create sub-mesh
//         fvMesh& cMesh = cMeshSubSet.subMesh();
//
//         // References for maps
//         const labelList& cCellMap = cMeshSubSet.cellMap();
//
//         // Create fields
//         volScalarField cSigma
//         (
//             IOobject
//             (
//                 "sigma",
//                 cMesh.time().timeName(),
//                 cMesh,
//                 IOobject::NO_READ
//             ),
//             cMesh,
//             dimensionedScalar(word(), sigma.dimensions(), pTraits<scalar>::zero),
//             zeroGradientFvPatchField<scalar>::typeName
//         );
//
//         volScalarField cAlpha
//         (
//             IOobject
//             (
//                 "alpha",
//                 cMesh.time().timeName(),
//                 cMesh,
//                 IOobject::NO_READ
//             ),
//             cMesh,
//             dimensionedScalar(word(), alpha.dimensions(), pTraits<scalar>::zero),
//             zeroGradientFvPatchField<scalar>::typeName
//         );
//
//         volVectorField cARe
//         (
//             IOobject
//             (
//                 "ARe",
//                 cMesh.time().timeName(),
//                 cMesh,
//                 IOobject::NO_READ
//             ),
//             cMesh,
//             dimensioned<vector>(word(), ARe.dimensions(), pTraits<vector>::zero),
//             zeroGradientFvPatchField<vector>::typeName
//         );
//
//         volVectorField cAIm
//         (
//             IOobject
//             (
//                 "AIm",
//                 cMesh.time().timeName(),
//                 cMesh,
//                 IOobject::NO_READ
//             ),
//             cMesh,
//             dimensionedVector(word(), AIm.dimensions(), vector::zero),
//             zeroGradientFvPatchField<vector>::typeName
//         );
//
//         volScalarField cVRe
//         (
//             IOobject
//             (
//                 "VRe",
//                 cMesh.time().timeName(),
//                 cMesh,
//                 IOobject::NO_READ
//             ),
//             cMesh,
//             dimensionedScalar(word(), VRe.dimensions(), pTraits<scalar>::zero),
//             zeroGradientFvPatchField<scalar>::typeName
//         );
//
//         volScalarField cVIm
//         (
//             IOobject
//             (
//                 "VIm",
//                 cMesh.time().timeName(),
//                 cMesh,
//                 IOobject::NO_READ
//             ),
//             cMesh,
//             dimensionedScalar(word(), VIm.dimensions(), pTraits<scalar>::zero),
//             zeroGradientFvPatchField<scalar>::typeName
//         );
//
//         // Map internal values of volume fields
//         // and correct their boundary values
//         {
//             forAll (cMesh.C(), cellI)
//             {
//                 cSigma.internalField()[cellI] = sigma.internalField()[cCellMap[cellI]];
//                 cAlpha.internalField()[cellI] = alpha.internalField()[cCellMap[cellI]];
//
//                 cARe.internalField()[cellI] = ARe.internalField()[cCellMap[cellI]];
//                 cAIm.internalField()[cellI] = AIm.internalField()[cCellMap[cellI]];
//
//                 cVRe.internalField()[cellI] = VRe.internalField()[cCellMap[cellI]];
//                 cVIm.internalField()[cellI] = VIm.internalField()[cCellMap[cellI]];
//             }
//
//             cSigma.correctBoundaryConditions();
//             cAlpha.correctBoundaryConditions();
//
//             cARe.correctBoundaryConditions();
//             cAIm.correctBoundaryConditions();
//
//             cVRe.correctBoundaryConditions();
//             cVIm.correctBoundaryConditions();
//         }
//
//         // Solve V equation
//         {
//             // Assemble V-equations
//             fvScalarMatrix cVReEqn
//             (
//                 fvm::laplacian(cSigma, cVRe, "laplacian(VRe)")
//               - fvc::div(cAlpha*cAIm, "div(AIm)")
//             );
//             fvScalarMatrix cVImEqn
//             (
//                 fvm::laplacian(cSigma, cVIm, "laplacian(VIm)")
//               + fvc::div(cAlpha*cARe, "div(ARe)")
//             );
//
//             label refCell = cMesh.findCell(point(0.00, 0.00, 0.01));
//             scalar refValue = 0.0;
//
//             cVReEqn.setReference(refCell, refValue);
//             cVImEqn.setReference(refCell, refValue);
//
//             // Solve V
//             cVReEqn.solve();
//             cVImEqn.solve();
//
//             cVRe.correctBoundaryConditions();
//             cVIm.correctBoundaryConditions();
//         }
//
//         // Reverse map fields
//         {
//             // Calculate gradient of electric scalar potential
//             volVectorField cVReGrad("VReGrad", fvc::grad(cVRe));
//             volVectorField cVImGrad("VImGrad", fvc::grad(cVIm));
//
//             // Calculate eddy currents
//             volVectorField cjRe("jRe",   cAlpha * cAIm - cSigma * cVReGrad);
//             volVectorField cjIm("jIm", - cAlpha * cARe - cSigma * cVImGrad);
//
//             forAll (cMesh.C(), cellI)
//             {
//                 VRe[cCellMap[cellI]] = cVRe[cellI];
//                 VIm[cCellMap[cellI]] = cVIm[cellI];
//
//                 VReGrad[cCellMap[cellI]] = cVReGrad[cellI];
//                 VImGrad[cCellMap[cellI]] = cVImGrad[cellI];
//
//                 jRe[cCellMap[cellI]] = cjRe[cellI];
//                 jIm[cCellMap[cellI]] = cjIm[cellI];
//             }
//         }
//     }

    jRe ==   alpha * AIm;
    jIm == - alpha * ARe;

    BRe == fvc::curl(ARe);
    BIm == fvc::curl(AIm);
