
    // Real part
    surfaceScalarField VReSource
    (
        "VReSource",
        sigmaf*fvc::interpolate(fvc::grad(VRe, "grad(VReSource)")) & mesh.Sf()
    );

    fvScalarMatrix VReEqn
    (
      - fvm::laplacian(sigmaf, VRe)
//       + fvc::div(alpha*AIm)
//      ==
//       - fvc::div(VReSource)
    );

    // Imaginary part
    surfaceScalarField VImSource
    (
        "VImSource",
        sigmaf*fvc::interpolate(fvc::grad(VIm, "grad(VImSource)")) & mesh.Sf()
    );

    fvScalarMatrix VImEqn
    (
      - fvm::laplacian(sigmaf, VIm)
//       - fvc::div(alpha*ARe)
//      ==
//       - fvc::div(VImSource)
    );

//     // Real part
//     fvScalarMatrix VReEqn
//     (
//       - fvm::laplacian(sigmaf, VRe)
// //       + fvc::div(alpha*AIm)
//     );
//
//     // Imaginary part
//     fvScalarMatrix VImEqn
//     (
//       - fvm::laplacian(sigmaf, VIm)
// //       - fvc::div(alpha*ARe)
//     );
