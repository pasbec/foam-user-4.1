
    // Assemble matrix for A
#   include "AEqn.H"

    // Assemble matrix for V
#   include "VEqn.H"

    // Solve AV-system
    if (AVmethod == "segregated")
    {
        int nSegregatedCorr = AVmethodDict.lookupOrDefault<int>("corr",10);

        for
        (
            int segregatedCorr=0;
            segregatedCorr<nSegregatedCorr;
            segregatedCorr++
        )
        {
#           include "solveAVsegregated.H"
        }
    }
    else if (AVmethod == "coupled")
    {
        int nCoupledCorr = AVmethodDict.lookupOrDefault<int>("corr",4);

        for
        (
            int coupledCorr=0;
            coupledCorr<nCoupledCorr;
            coupledCorr++
        )
        {
#           include "solveAVcoupled.H"
        }
    }
    else
    {
        FatalErrorIn("AVEqn.H")
        << "Method for solving AV-system '" << AVmethod
        << "' is not supported"
        << abort(FatalError);
    }

    // Update corrected gradient of V
    {
        // Interpolate V to face centers
        surfaceScalarField sVRe
        (
            "sVRe",
            fvc::interpolate(VRe)
        );
        surfaceScalarField sVIm
        (
            "sVIm",
            fvc::interpolate(VIm)
        );

        // Enforce zero normal gradient for
        // V (dV/dn=0)
        {
            // Mesh cell-face addressing
            const unallocLabelList& own = mesh.owner();
            const unallocLabelList& nei = mesh.neighbour();

            // Replace face values with cell values
            // from conducting side
            forAll (cnPNfaces, facei)
            {
                label faceI = cnPNfaces[facei];
                sVRe[faceI] = VRe[own[faceI]];
                sVIm[faceI] = VIm[own[faceI]];
            }
            forAll (cnNPfaces, facei)
            {
                label faceI = cnNPfaces[facei];
                sVRe[faceI] = VRe[nei[faceI]];
                sVIm[faceI] = VIm[nei[faceI]];
            }
        }

        // Update corrected gradient
        VReGrad = fvc::grad(sVRe);
        VImGrad = fvc::grad(sVIm);

        // Remove all gradient values from
        // non-conducting region
        VReGrad *= sigmaIc;
        VImGrad *= sigmaIc;
    }



    // Calculate derived fields
    {
        // Magnetic field density
        BRe == fvc::curl(ARe);
        BIm == fvc::curl(AIm);

        // Eddy current density
        jRe ==   alpha * AIm - sigma * VReGrad;
        jIm == - alpha * ARe - sigma * VImGrad;
    }



//     // Gauging of A
//     if (!gauged)
//     {
//         AReDiv == fvc::div(ARe);
//         AImDiv == fvc::div(AIm);
//
//         gauged = true;
//     }

//     // Gauging of A
//     if (!gauged)
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
//         gauged = true;
//     }




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

//     // Gauge fixing: Solve for gauge fixing potential G
//     {
//         // Assemble G-equations
//         fvScalarMatrix GReEqn
//         (
//             fvm::laplacian(GRe)
//           + AReDiv
//         );
//         fvScalarMatrix GImEqn
//         (
//             fvm::laplacian(GIm)
//           + AImDiv
//         );
//
//         // Reference settings
//         label refCell = mesh.findCell(point(0.00, 0.00, 0.01));
//         scalar refValue = 0.0;
//
//         GReEqn.setReference(refCell, refValue);
//         GImEqn.setReference(refCell, refValue);
//
//         // Solve G
//         GReEqn.solve();
//         GImEqn.solve();
//
//         GRe.correctBoundaryConditions();
//         GIm.correctBoundaryConditions();
//
//         // TODO: Remove after debug
//         volVectorField AReCoulomb("AReCoulomb", ARe + fvc::grad(GRe));
//         volVectorField AImCoulomb("AImCoulomb", AIm + fvc::grad(GIm));
//         AReCoulomb.write();
//         AImCoulomb.write();
//     }
//     // Gauge fixing: Solve for gauge fixing potential G
//     {
//         // Assemble G-equations
//         fvScalarMatrix GReEqn
//         (
//             fvm::laplacian(GRe)
//           + fvc::div(ARe)
//         );
//         fvScalarMatrix GImEqn
//         (
//             fvm::laplacian(GIm)
//           + fvc::div(AIm)
//         );
//
//         // Modify fvm::laplacian() for GRe and GIm to enforce
//         // zero normal gradient (dGRe/dn=0) and (dGIm/dn=0) at
//         // the conductor interface as an inner constraint:
//         //
//         // dG/dn=0: G_N - G_P / |d| = G_f = 0
//         //
//         // dG/dn=0: G_P = G_N
//         //          G_N = G_P
//         //
//         // dG/dn=0: diag_P = diag_P + upper
//         //          diag_N = diag_N + lower
//         //          upper = 0
//         //          lower = 0
//         //
//         {
//             // Mesh cell-face addressing
//             const unallocLabelList& own = mesh.owner();
//             const unallocLabelList& nei = mesh.neighbour();
//
//             // Loop over all faces of conductor interface
//             forAll (cnFaces, facei)
//             {
//                 label faceI = cnFaces[facei];
//
//                 // Eliminate upper/lower coeffs based on constraint
//                 // NOTE: This will keep SYMMETRIC matrices and it will
//                 //       also keep their diagonal dominance
//                 scalar GReEqnUf = GReEqn.upper()[faceI];
//                 GReEqn.diag()[own[faceI]] += GReEqnUf;
//                 GReEqn.diag()[nei[faceI]] += GReEqnUf;
//                 GReEqn.upper()[faceI] -= GReEqnUf;
//                 scalar GImEqnUf = GImEqn.upper()[faceI];
//                 GImEqn.diag()[own[faceI]] += GImEqnUf;
//                 GImEqn.diag()[nei[faceI]] += GImEqnUf;
//                 GImEqn.upper()[faceI] -= GImEqnUf;
//
//                 // TODO: What about gradient correction?
//                 // We have to modifiy the source term, too. But
//                 // this depends on the type of correction. Hmm...
//             }
//         }
//
//         // Reference settings
//         label refCell = mesh.findCell(point(0.00, 0.00, 0.01));
//         scalar refValue = 0.0;
//
//         GReEqn.setReference(refCell, refValue);
//         GImEqn.setReference(refCell, refValue);
//
//         // Solve G
//         GReEqn.solve();
//         GImEqn.solve();
//
//         GRe.correctBoundaryConditions();
//         GIm.correctBoundaryConditions();
//
//         // Interpolate G to face centers
//         surfaceScalarField csGRe
//         (
//             "csGRe",
//             fvc::interpolate(GRe)
//         );
//         surfaceScalarField csGIm
//         (
//             "csGIm",
//             fvc::interpolate(GIm)
//         );
//         surfaceScalarField nsGRe
//         (
//             "nsGRe",
//             fvc::interpolate(GRe)
//         );
//         surfaceScalarField nsGIm
//         (
//             "nsGIm",
//             fvc::interpolate(GIm)
//         );
//
//         // Update corrected gradient of G
//         {
//             // Enforce zero normal gradient for
//             // G (dG/dn=0)
//             {
//                 // Mesh cell-face addressing
//                 const unallocLabelList& own = mesh.owner();
//                 const unallocLabelList& nei = mesh.neighbour();
//
//                 // Replace face values with cell values
//                 // from conducting/non-conducting side
//                 forAll (cnPNfaces, facei)
//                 {
//                     label faceI = cnPNfaces[facei];
//                     csGRe[faceI] = GRe[own[faceI]];
//                     csGIm[faceI] = GIm[own[faceI]];
//                     nsGRe[faceI] = GRe[nei[faceI]];
//                     nsGIm[faceI] = GIm[nei[faceI]];
//                 }
//                 forAll (cnNPfaces, facei)
//                 {
//                     label faceI = cnNPfaces[facei];
//                     csGRe[faceI] = GRe[nei[faceI]];
//                     csGIm[faceI] = GIm[nei[faceI]];
//                     nsGRe[faceI] = GRe[own[faceI]];
//                     nsGIm[faceI] = GIm[own[faceI]];
//                 }
//             }
//
//             // Update corrected gradient
//             GReGrad = sigmaIc*fvc::grad(csGRe) + sigmaIn*fvc::grad(nsGRe);
//             GImGrad = sigmaIc*fvc::grad(csGIm) + sigmaIn*fvc::grad(nsGIm);
//         }
//
//         // TODO: Remove after debug
//         volVectorField AReCoulomb("AReCoulomb", ARe + GReGrad);
//         volVectorField AImCoulomb("AImCoulomb", AIm + GImGrad);
//         AReCoulomb.write();
//         AImCoulomb.write();
//     }
