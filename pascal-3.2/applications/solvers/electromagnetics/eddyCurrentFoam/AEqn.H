
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
