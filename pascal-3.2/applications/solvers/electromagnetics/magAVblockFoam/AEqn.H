
    // Real part
    fvVectorMatrix AReEqn
    (
      - rMu0 * fvm::laplacian(ARe, "laplacian(ARe)")
//       - alpha * AIm
//       + sigma * fvc::grad(VRe)
     ==
        jsRe
    );



    // Imaginary part
    fvVectorMatrix AImEqn
    (
      - rMu0 * fvm::laplacian(AIm, "laplacian(AIm)")
//       + alpha * ARe
//       + sigma * fvc::grad(VIm)
     ==
        jsIm
    );
