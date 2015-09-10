
    // Real part
    fvVectorMatrix AReEqn
    (
        scaleA * rMu0 * fvm::laplacian(ARe)
//       + scaleA * alpha * AIm
//       - scaleA * sigma * fvc::grad(VRe)
      + scaleA * jsRe
    );

    // Imaginary part
    fvVectorMatrix AImEqn
    (
        scaleA * rMu0 * fvm::laplacian(AIm)
//       - scaleA * alpha * ARe
//       - scaleA * sigma * fvc::grad(VIm)
      + scaleA * jsIm
    );
