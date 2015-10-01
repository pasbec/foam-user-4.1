
    // Create list with boundary types
    wordList ARePatches(mesh.boundary().size());
    wordList AImPatches(mesh.boundary().size());

    forAll(ARe.boundaryField(), patchI)
    {
        ARePatches[patchI] =
            ARe.boundaryField()[patchI].type();
    }
    forAll(AIm.boundaryField(), patchI)
    {
        AImPatches[patchI] =
            AIm.boundaryField()[patchI].type();
    }

    // TODO: Check, that ARe/AIm use correct boundary
    //       conditions! Currently, only tangentialMagnetic
    //       is supported!

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
        ARe.component(0),
        ARePatches
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
        ARe.component(1),
        ARePatches
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
        ARe.component(2),
        ARePatches
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
        AIm.component(0),
        AImPatches
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
        AIm.component(1),
        AImPatches
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
        AIm.component(2),
        AImPatches
    );



    // Assemble A-equations
    fvScalarMatrix ARexEqn
    (
        rMu0 * fvm::laplacian(ARex)
//       - rMu0 * AReDivGrad.component(0)
//    + alpha * AImx
//    - sigma * VReGrad.component(0)
      + jsRe.component(0)
    );
    fvScalarMatrix AReyEqn
    (
        rMu0 * fvm::laplacian(ARey)
//       - rMu0 * AReDivGrad.component(1)
//    + alpha * AImy
//    - sigma * VReGrad.component(1)
      + jsRe.component(1)
    );
    fvScalarMatrix ARezEqn
    (
        rMu0 * fvm::laplacian(ARez)
//       - rMu0 * AReDivGrad.component(2)
//    + alpha * AImz
//    - sigma * VReGrad.component(2)
      + jsRe.component(2)
    );
    fvScalarMatrix AImxEqn
    (
        rMu0 * fvm::laplacian(AImx)
//       - rMu0 * AImDivGrad.component(0)
//    - alpha * ARex
//    - sigma * VImGrad.component(0)
      + jsIm.component(0)
    );
    fvScalarMatrix AImyEqn
    (
        rMu0 * fvm::laplacian(AImy)
//       - rMu0 * AImDivGrad.component(1)
//    - alpha * ARey
//    - sigma * VImGrad.component(1)
      + jsIm.component(1)
    );
    fvScalarMatrix AImzEqn
    (
        rMu0 * fvm::laplacian(AImz)
//       - rMu0 * AImDivGrad.component(2)
//    - alpha * ARez
//    - sigma * VImGrad.component(2)
      + jsIm.component(2)
    );



//     // Modify fvm::laplacian() for all components of ARe and AIm to
//     // enforce zero normal components (ARe*n=0) and (AIm*n=0) at the
//     // conductor interface as an inner constraint:
//     //
//     // A*Sf=Af=0: Af_P = - w/(1-w) * Af_N
//     //            Af_N = - (1-w)/w * Af_P
//     //
//     // A*Sf=Af=0: diag_P = TODO
//     //            diag_N = TODO
//     //            upper = 0
//     //            lower = 0
//     {
//         // Mesh cell-face addressing
//         const unallocLabelList& own = mesh.owner();
//         const unallocLabelList& nei = mesh.neighbour();
//
//         // Mesh geometry
//         const vectorField& Cf = mesh.Cf();
//         const vectorField& C = mesh.C();
//         const vectorField& Sf = mesh.Sf();
//         const scalarField& magSf = mesh.magSf();
//
//         // Derived mesh geometry
//         const vectorField Sn = Sf/magSf;
//         const vectorField SnCmptMag = cmptMag(Sn);
//         const scalarField SnxMag = mag(Sn.component(0));
//         const scalarField SnyMag = mag(Sn.component(1));
//         const scalarField SnzMag = mag(Sn.component(2));
//
//         // Loop over all faces of conductor interface
//         forAll (cnFaces, facei)
//         {
//             label faceI = cnFaces[facei];
//
//             // Calculate face centre/cell centre distances Pf and fN
//             scalar Pf = mag(Sn[faceI] & (Cf[faceI] - C[own[faceI]]));
//             scalar fN = mag(Sn[faceI] & (C[nei[faceI]] - Cf[faceI]));
//
//             // Calculate coefficient constraint weights
//             scalar w = fN / (Pf+fN);
//             scalar wOwn = w/(1-w);
//             scalar wNei = (1-w)/w;
//
//             // Calculate coefficient direction weights
//             scalar wSnx = SnxMag[faceI];
//             scalar wSny = SnyMag[faceI];
//             scalar wSnz = SnzMag[faceI];
//
//             // Eliminate upper/lower coeffs based on constraint
//             // NOTE: This will keep SYMMETRIC matrices and it will
//             //       also sligthly improve their diagonal dominance
//             scalar wSnARexEqnUf = wSnx * ARexEqn.upper()[faceI];
//             ARexEqn.diag()[own[faceI]] -= wSnARexEqnUf * wOwn;
//             ARexEqn.diag()[nei[faceI]] -= wSnARexEqnUf * wNei;
//             ARexEqn.upper()[faceI] -= wSnARexEqnUf;
//             scalar wSnAReyEqnUf = wSny * AReyEqn.upper()[faceI];
//             AReyEqn.diag()[own[faceI]] -= wSnAReyEqnUf * wOwn;
//             AReyEqn.diag()[nei[faceI]] -= wSnAReyEqnUf * wNei;
//             AReyEqn.upper()[faceI] -= wSnAReyEqnUf;
//             scalar wSnARezEqnUf = wSnz * ARezEqn.upper()[faceI];
//             ARezEqn.diag()[own[faceI]] -= wSnARezEqnUf * wOwn;
//             ARezEqn.diag()[nei[faceI]] -= wSnARezEqnUf * wNei;
//             ARezEqn.upper()[faceI] -= wSnARezEqnUf;
//             scalar wSnAImxEqnUf = wSnx * AImxEqn.upper()[faceI];
//             AImxEqn.diag()[own[faceI]] -= wSnAImxEqnUf * wOwn;
//             AImxEqn.diag()[nei[faceI]] -= wSnAImxEqnUf * wNei;
//             AImxEqn.upper()[faceI] -= wSnAImxEqnUf;
//             scalar wSnAImyEqnUf = wSny * AImyEqn.upper()[faceI];
//             AImyEqn.diag()[own[faceI]] -= wSnAImyEqnUf * wOwn;
//             AImyEqn.diag()[nei[faceI]] -= wSnAImyEqnUf * wNei;
//             AImyEqn.upper()[faceI] -= wSnAImyEqnUf;
//             scalar wSnAImzEqnUf = wSnz * AImzEqn.upper()[faceI];
//             AImzEqn.diag()[own[faceI]] -= wSnAImzEqnUf * wOwn;
//             AImzEqn.diag()[nei[faceI]] -= wSnAImzEqnUf * wNei;
//             AImzEqn.upper()[faceI] -= wSnAImzEqnUf;
//
//             // TODO: What about gradient correction?
//             // We have to modifiy the source term, too. But
//             // this depends on the type of correction. Hmm...
//         }
//     }



//     // Add explicit terms from gauge G
//     {
//         volScalarField GReGradDiv
//         (
//             "GReGradDiv",
//             fvc::div(GReGrad)
//         );
//         volScalarField GImGradDiv
//         (
//             "GImGradDiv",
//             fvc::div(GImGrad)
//         );
//
//         volVectorField GReGradDivGrad
//         (
//             "GReGradDivGrad",
//             fvc::grad(AReDiv)
//         );
//         volVectorField GImGradDivGrad
//         (
//             "GImGradDivGrad",
//             fvc::grad(AImDiv)
//         );
//
//         ARexEqn += rMu0 * GReGradDivGrad.component(0);
//         AReyEqn += rMu0 * GReGradDivGrad.component(1);
//         ARezEqn += rMu0 * GReGradDivGrad.component(2);
//         AImxEqn += rMu0 * GImGradDivGrad.component(0);
//         AImyEqn += rMu0 * GImGradDivGrad.component(1);
//         AImzEqn += rMu0 * GImGradDivGrad.component(2);
//
//         // TODO: Remove after debug
//         GReGradDiv.write();
//         GImGradDiv.write();
//         GReGradDivGrad.write();
//         GImGradDivGrad.write();
//     }



//     // Add explicit terms which are necessary if div(A) is NOT zero
//     {
//         // Interpolate A to face centers
//         surfaceScalarField sARe
//         (
//             "sARe",
//             fvc::interpolate(ARe) & mesh.Sf()
//         );
//         surfaceScalarField sAIm
//         (
//             "sAIm",
//             fvc::interpolate(AIm) & mesh.Sf()
//         );
//
//         // Enforce enforce zero normal components (ARe*n=0)
//         // and (AIm*n=0) at the conductor interface as an
//         // inner constraint
//         {
//             // Replace face values with zeros
//             forAll (cnFaces, facei)
//             {
//                 label faceI = cnFaces[facei];
//                 sARe[faceI] = 0.0;
//                 sAIm[faceI] = 0.0;
//             }
//         }
//
//         AReDiv == fvc::div(sARe);
//         AImDiv == fvc::div(sAIm);
//
// //         AReDiv == fvc::div(ARe);
// //         AImDiv == fvc::div(AIm);
//
// //         volScalarField AReDiv
// //         (
// //             "AReDiv",
// //             fvc::div(ARe)
// //         );
// //         volScalarField AImDiv
// //         (
// //             "AImDiv",
// //             fvc::div(AIm)
// //         );
//
//         volVectorField AReDivGrad
//         (
//             "AReDivGrad",
//             fvc::grad(AReDiv)
//         );
//         volVectorField AImDivGrad
//         (
//             "AImDivGrad",
//             fvc::grad(AImDiv)
//         );
//
//         ARexEqn -= rMu0 * AReDivGrad.component(0);
//         AReyEqn -= rMu0 * AReDivGrad.component(1);
//         ARezEqn -= rMu0 * AReDivGrad.component(2);
//         AImxEqn -= rMu0 * AImDivGrad.component(0);
//         AImyEqn -= rMu0 * AImDivGrad.component(1);
//         AImzEqn -= rMu0 * AImDivGrad.component(2);
//
//         // TODO: Remove after debug
//         AReDiv.write();
//         AImDiv.write();
//         AReDivGrad.write();
//         AImDivGrad.write();
//     }
