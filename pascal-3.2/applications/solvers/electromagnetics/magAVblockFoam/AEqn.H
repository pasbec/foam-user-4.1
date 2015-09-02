
    // Real part
    fvVectorMatrix AReEqn
    (
      - rMu0 * fvm::laplacian(ARe, "laplacian(ARe)")
//       + rMu0 * fvc::grad(divARe)
//       - alpha * AIm
//       + sigma * fvc::grad(VRe)
     ==
        jsRe
    );

//     AReEqn.relax(0.9999);



    // Imaginary part
    fvVectorMatrix AImEqn
    (
      - rMu0 * fvm::laplacian(AIm, "laplacian(AIm)")
//       + rMu0 * fvc::grad(divAIm)
//       + alpha * ARe
//       + sigma * fvc::grad(VIm)
     ==
        jsIm
    );

//     AImEqn.relax(0.9999);
