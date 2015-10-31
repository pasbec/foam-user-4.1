
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
        scalarField& VImEqnD = VImEqn.diag();
        scalarField& VImEqnU = VImEqn.upper();

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
    }
