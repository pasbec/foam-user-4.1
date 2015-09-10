
    // Real part
    fvVectorMatrix AReEqn
    (
        rMu0 * fvm::laplacian(ARe)
//       + alpha * AIm
//       - sigma * fvc::grad(VRe)
      + jsRe
    );

    // Imaginary part
    fvVectorMatrix AImEqn
    (
        rMu0 * fvm::laplacian(AIm)
//       - alpha * ARe
//       - sigma * fvc::grad(VIm)
      + jsIm
    );
