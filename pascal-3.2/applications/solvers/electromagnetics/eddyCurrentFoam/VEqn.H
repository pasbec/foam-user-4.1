
    // Get reference to sigmaV pointer
    volScalarField& sigmaV = *sigmaVPtr;

    // Interpolate sigma to face centers
    surfaceScalarField sSigmaV
    (
        "sSigmaV",
        fvc::interpolate(sigmaV)
    );

    // Enforce zero normal gradient for sigmaV
    // (dsigmaV/dn=0) at the conductor interface as
    // inner constrain
    {
        // Mesh cell-face addressing
        const unallocLabelList& own = mesh.owner();
        const unallocLabelList& nei = mesh.neighbour();

        // Replace face values with cell values
        // from conducting side
        forAll (cnPNfaces, facei)
        {
            label faceI = cnPNfaces[facei];
            sSigmaV[faceI] = sigmaV[own[faceI]];
        }
        forAll (cnNPfaces, facei)
        {
            label faceI = cnNPfaces[facei];
            sSigmaV[faceI] = sigmaV[nei[faceI]];
        }
    }

    // Calculate gradient of alpha
    volVectorField gradAlpha
    (
        "gradAlpha",
        omega * fvc::grad(sSigmaV)
    );

    // Assemble V-equations
    fvScalarMatrix VReEqn
    (
        fvm::laplacian(sSigmaV, VRe, "laplacian(sigmaV,VRe)")
//    - fvc::div(sigmaV * AIm)
    );
    fvScalarMatrix VImEqn
    (
        fvm::laplacian(sSigmaV, VIm, "laplacian(sigmaV,VIm)")
//    + fvc::div(sigmaV * ARe)
    );

    fvScalarMatrix VReEqnCorr = VReEqn;
    fvScalarMatrix VImEqnCorr = VImEqn;

    // Temporarily force zero gradient for V (dV/dn=0)
    // at conductor interface as inner constraint. We
    // need to do this in order set a new gradient
    // for V depnding on A, later!
    {
        // Grab parts of V-equations
        scalarField& VReEqnD = VReEqn.diag();
        scalarField& VReEqnU = VReEqn.upper();
        scalarField& VReEqnS = VReEqn.source();
        scalarField& VImEqnD = VImEqn.diag();
        scalarField& VImEqnU = VImEqn.upper();
        scalarField& VImEqnS = VImEqn.source();

        // Mesh cell-face addressing
        const unallocLabelList& own = mesh.owner();
        const unallocLabelList& nei = mesh.neighbour();

        // Loop over all faces of conductor interface
        forAll (cnFaces, facei)
        {
            label faceI = cnFaces[facei];

            // Write
            {
                VReEqnD[own[faceI]] += VReEqnU[faceI];
                VReEqnD[nei[faceI]] += VReEqnU[faceI];
                VReEqnU[faceI]      -= VReEqnU[faceI];

                VImEqnD[own[faceI]] += VImEqnU[faceI];
                VImEqnD[nei[faceI]] += VImEqnU[faceI];
                VImEqnU[faceI]      -= VImEqnU[faceI];
            }

            // TODO: What about gradient correction?
            // We have to modifiy the source term, too. But
            // this depends on the type of correction. Hmm...
        }

        // Reference settings
        {
            // Reference values in conducting region
            scalar cRefValueVRe = 0.0;
            scalar cRefValueVIm = 0.0;

            // Force reference cells in conducting region
            VReEqnD[cRefCell] += VReEqnD[cRefCell];
            VImEqnD[cRefCell] += VImEqnD[cRefCell];
            VReEqnS[cRefCell] += VReEqnD[cRefCell]*cRefValueVRe;
            VImEqnS[cRefCell] += VImEqnD[cRefCell]*cRefValueVIm;

//             // TEST
//             // FIXME: Refpoint problems!
//             // NOTE: Pseudo time-stepping based on Euler-explicit
//             // NOTE: Minus sign as laplacian is in matrix with plus sign
//             scalarField V = mesh.V();
//             scalar pseudoDeltaT = 1.0;
//             VReEqnD -= V/pseudoDeltaT;
//             VImEqnD -= V/pseudoDeltaT;
//             VReEqnS -= V/pseudoDeltaT * (VRe.internalField() - average(VRe.internalField()));
//             VImEqnS -= V/pseudoDeltaT * (VIm.internalField() - average(VIm.internalField()));
//             VReEqnD += 1.0e-04 * VReEqnD;
//             VImEqnD += 1.0e-04 * VImEqnD;
//             VReEqnS += 1.0e-04 * VReEqnD * (VRe.internalField() - average(VRe.internalField()));
//             VImEqnS += 1.0e-04 * VImEqnD * (VIm.internalField() - average(VIm.internalField()));

//             // TEST
//             // FIXME: Refpoint problems!
//             VReEqnD[cRefCell] += 0.01*VReEqnD[cRefCell];
//             VImEqnD[cRefCell] += 0.01*VImEqnD[cRefCell];
//             VReEqnS[cRefCell] += 0.01*VReEqnD[cRefCell]*cRefValueVRe;
//             VImEqnS[cRefCell] += 0.01*VImEqnD[cRefCell]*cRefValueVIm;

//             // TEST
//             // FIXME: Refpoint problems!
//             VReEqn.setValues(labelList(1,cRefCell), scalarField(1,cRefValueVRe));
//             VImEqn.setValues(labelList(1,cRefCell), scalarField(1,cRefValueVIm));

//             // TEST
//             // FIXME: Refpoint problems!
//             Info << "max(VReEqn.diag()) = " << max(VReEqn.diag()) << endl;
//             Info << "min(VReEqn.diag()) = " << min(VReEqn.diag()) << endl;
//             Info << "max(VReEqn.source()) = " << max(VReEqn.source()) << endl;
//             Info << "min(VReEqn.source()) = " << min(VReEqn.source()) << endl;
//
//             Info << "VReEqn.diag()[cRefCell] = " << VReEqn.diag()[cRefCell] << endl;
//             Info << "VReEqn.diag()[cRefCell] = " << VReEqn.diag()[cRefCell] << endl;
//             Info << "VReEqn.source()[cRefCell] = " << VReEqn.source()[cRefCell] << endl;
//             Info << "VReEqn.source()[cRefCell] = " << VReEqn.source()[cRefCell] << endl;
//
//             Info << "max(VImEqn.diag()) = " << max(VImEqn.diag()) << endl;
//             Info << "min(VImEqn.diag()) = " << min(VImEqn.diag()) << endl;
//             Info << "max(VImEqn.source()) = " << max(VImEqn.source()) << endl;
//             Info << "min(VImEqn.source()) = " << min(VImEqn.source()) << endl;
//
//             Info << "VImEqn.diag()[cRefCell] = " << VImEqn.diag()[cRefCell] << endl;
//             Info << "VImEqn.diag()[cRefCell] = " << VImEqn.diag()[cRefCell] << endl;
//             Info << "VImEqn.source()[cRefCell] = " << VImEqn.source()[cRefCell] << endl;
//             Info << "VImEqn.source()[cRefCell] = " << VImEqn.source()[cRefCell] << endl;

            // Reference values in non-conducting region
            // TODO: Which values to use here?
            scalar nRefValueVRe = 0.0;
            scalar nRefValueVIm = 0.0;

            // Set reference cells in non-conducting region if needed
            VReEqn.setReference(nRefCell, nRefValueVRe);
            VImEqn.setReference(nRefCell, nRefValueVIm);
        }
    }


//     // TEST
//     // FIXME: Refpoint problems!
//     {
//         // Grab parts of V-equations
//         scalarField& VReEqnD = VReEqnCorr.diag();
//         scalarField& VReEqnU = VReEqnCorr.upper();
//         scalarField& VReEqnS = VReEqnCorr.source();
//         scalarField& VImEqnD = VImEqnCorr.diag();
//         scalarField& VImEqnU = VImEqnCorr.upper();
//         scalarField& VImEqnS = VImEqnCorr.source();
//
//         // Mesh cell-face addressing
//         const unallocLabelList& own = mesh.owner();
//         const unallocLabelList& nei = mesh.neighbour();
//
//         // Mesh geometry
//         const vectorField& Sf = mesh.Sf();
//         const scalarField& w = mesh.weights();
//
//         // Loop over all faces of conductor interface
//         forAll (cnFaces, facei)
//         {
//             label faceI = cnFaces[facei];
//
//             // Calculate weight coeffs
//             scalar wP = omega.value() * w[faceI];
//             scalar wN = omega.value() - wP;
//
//             // Write
//             {
//                 VReEqnD[own[faceI]] += VReEqnU[faceI];
//                 VReEqnD[nei[faceI]] += VReEqnU[faceI];
//                 VReEqnU[faceI]      -= VReEqnU[faceI];
//
//                 VImEqnD[own[faceI]] += VImEqnU[faceI];
//                 VImEqnD[nei[faceI]] += VImEqnU[faceI];
//                 VImEqnU[faceI]      -= VImEqnU[faceI];
//             }
//
//             // TODO: What about gradient correction?
//             // We have to modifiy the source term, too. But
//             // this depends on the type of correction. Hmm...
//         }
//     }
