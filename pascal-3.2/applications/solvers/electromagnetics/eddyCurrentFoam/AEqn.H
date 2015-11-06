
    // Assemble base-matrices for A
    fvVectorMatrix AReEqn
    (
        rMu0 * fvm::laplacian(ARe)
//    + alpha * AIm
//    - sigma * VReGrad
      + jsRe
    );
    fvVectorMatrix AImEqn
    (
        rMu0 * fvm::laplacian(AIm)
//    - alpha * ARe
//    - sigma * VImGrad
      + jsIm
    );



    // Create list with boundary types
    wordList ARePatchTypes(mesh.boundary().size());
    wordList AImPatchTypes(mesh.boundary().size());

    forAll(ARe.boundaryField(), patchI)
    {
        ARePatchTypes[patchI] =
            ARe.boundaryField()[patchI].type();
    }
    forAll(AIm.boundaryField(), patchI)
    {
        AImPatchTypes[patchI] =
            AIm.boundaryField()[patchI].type();
    }
