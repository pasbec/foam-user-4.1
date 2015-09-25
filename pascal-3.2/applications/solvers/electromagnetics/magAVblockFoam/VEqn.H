
    // Get reference to sigmaV/sigmaI pointer
    volScalarField& sigmaV = *sigmaVPtr;
    volScalarField& sigmaIc = *sigmaIcPtr;
    volScalarField& sigmaIn = *sigmaInPtr;

    volScalarField rSigmaV
    (
        "rSigmaV",
        1.0/sigmaV
    );

    // Interpolate sigma to face centers
    surfaceScalarField sSigmaV
    (
        "sSigmaV",
        fvc::interpolate(sigmaV)
    );

    // Enforce a zero normal gradient for
    // sigmaV (d(sigmaV)/dn=0)
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

    // Create scale factor for V-equation in relation
    // to A-equation based on dimensional relation
    // [A] <=> [grad(V)]
    dimensionedScalar scaleV
    (
        "scaleV",
        pow(average(mesh.V()), scalar(1.0/3.0))
    );

    // Assemble V-equations
    fvScalarMatrix VReEqn
    (
        scaleV * fvm::laplacian(sSigmaV, VRe, "laplacian(sigmaV,VRe)")
//    - scaleV * omega * fvc::div(sigmaV * AIm)
    );
    fvScalarMatrix VImEqn
    (
        scaleV * fvm::laplacian(sSigmaV, VIm, "laplacian(sigmaV,VIm)")
//    + scaleV * omega * fvc::div(sigmaV * ARe)
    );

    // Modify fvm::laplacian() for VRe and VIm to enforce
    // zero normal gradient (dVRe/dn=0) and (dVIm/dn=0) at
    // the conductor interface as an inner constraint
    {
        // Mesh cell-face addressing
        const unallocLabelList& own = mesh.owner();
        const unallocLabelList& nei = mesh.neighbour();

        // Loop over all faces of conductor interface
        forAll (cnFaces, facei)
        {
            label faceI = cnFaces[facei];

            // Eliminate upper/lower coeffs based on constraint
            // NOTE: This will keep SYMMETRIC matrices and it will
            //       also keep their diagonal dominance
            scalar VReEqnUf = VReEqn.upper()[faceI];
            VReEqn.diag()[own[faceI]] += VReEqnUf;
            VReEqn.diag()[nei[faceI]] += VReEqnUf;
            VReEqn.upper()[faceI] -= VReEqnUf;
            scalar VImEqnUf = VImEqn.upper()[faceI];
            VImEqn.diag()[own[faceI]] += VImEqnUf;
            VImEqn.diag()[nei[faceI]] += VImEqnUf;
            VImEqn.upper()[faceI] -= VImEqnUf;

            // TODO: What about gradient correction?
            // We have to modifiy the source term, too. But
            // this depends on the type of correction. Hmm...
        }
    }
