
    // Real part
    fvVectorMatrix AReEqn
    (
      - rMu0 * fvm::laplacian(ARe)
//       + rMu0 * fvc::grad(divARe)
//       - alpha * AIm
//       + sigma * fvc::grad(VRe)
     ==
        jsRe
    );

    // Imaginary part
    fvVectorMatrix AImEqn
    (
      - rMu0 * fvm::laplacian(AIm)
//       + rMu0 * fvc::grad(divAIm)
//       + alpha * ARe
//       + sigma * fvc::grad(VIm)
     ==
        jsIm
    );
