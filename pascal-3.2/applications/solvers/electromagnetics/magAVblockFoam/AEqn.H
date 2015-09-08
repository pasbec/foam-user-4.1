
    // Real part
    fvVectorMatrix AReEqn
    (
        Ascale * rMu0 * fvm::laplacian(ARe, "laplacian(ARe)")
//       + Ascale * alpha * AIm
//       - Ascale * sigma * fvc::grad(VRe)
     ==
      - Ascale * jsRe
    );



    // Imaginary part
    fvVectorMatrix AImEqn
    (
        Ascale * rMu0 * fvm::laplacian(AIm, "laplacian(AIm)")
//       - Ascale * alpha * ARe
//       - Ascale * sigma * fvc::grad(VIm)
     ==
      - Ascale * jsIm
    );
