
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

    // Temporarily force zero gradient for V (dV/dn=0)
    // at conductor interface as inner constraint. We
    // need to do this in order set a new gradient
    // for V depnding on A, later!
    {
        // Mesh cell-face addressing
        const unallocLabelList& own = mesh.owner();
        const unallocLabelList& nei = mesh.neighbour();

        // Loop over all faces of conductor interface
        forAll (cnFaces, facei)
        {
            label faceI = cnFaces[facei];

            // Write
            {
                scalar VReEqnUf = VReEqn.upper()[faceI];
                VReEqn.upper()[faceI] -= VReEqnUf;
                VReEqn.diag()[own[faceI]] += VReEqnUf;
                VReEqn.diag()[nei[faceI]] += VReEqnUf;

                scalar VImEqnUf = VImEqn.upper()[faceI];
                VImEqn.upper()[faceI] -= VImEqnUf;
                VImEqn.diag()[own[faceI]] += VImEqnUf;
                VImEqn.diag()[nei[faceI]] += VImEqnUf;
            }

            // TODO: What about gradient correction?
            // We have to modifiy the source term, too. But
            // this depends on the type of correction. Hmm...
        }

        // Reference settings
        {
            // Reference face
            label refFace = cnPNfaces[0];

            // Reference values in conducting region
            scalar refValueCVRe = 0;
            scalar refValueCVIm = 0;

            // Force reference cell in conducting region (owner)
            VReEqn.diag()[own[refFace]] += VReEqn.diag()[own[refFace]];
            VImEqn.diag()[own[refFace]] += VImEqn.diag()[own[refFace]];
            VReEqn.source()[own[refFace]] += VReEqn.diag()[own[refFace]]*refValueCVRe;
            VImEqn.source()[own[refFace]] += VImEqn.diag()[own[refFace]]*refValueCVIm;

            // Reference values in non-conducting region (neighbour)
            // TODO: Which values to use here?
            scalar refValueNVRe = 0;
            scalar refValueNVIm = 0;

            // Set reference cell in non-conducting region if needed
            VReEqn.setReference(nei[refFace], refValueNVRe);
            VImEqn.setReference(nei[refFace], refValueNVIm);
        }
    }
