
    // Real part
    fvScalarMatrix GReEqn
    (
      - fvm::laplacian(GRe, "laplacian(GRe)")
//       - fvc::div(ARe)
    );



    // Imaginary part
    fvScalarMatrix GImEqn
    (
      - fvm::laplacian(GIm, "laplacian(GIm)")
//       - fvc::div(AIm)
    );
