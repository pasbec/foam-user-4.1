
//     // Real part
//     surfaceScalarField rAReAf
//     (
//         "rAReAf",
//         sigmaf*fvc::interpolate(1.0/AReEqn.A())
//     );
// 
//     surfaceScalarField VReSource
//     (
//         "VReSource",
//         rAReAf*(fvc::interpolate(fvc::grad(VRe, "grad(VReSource)")) & mesh.Sf())
//     );
//     
//     fvScalarMatrix VReEqn
//     (
//       - fvm::laplacian(rAReAf, VRe)
//      ==
// //       - fvc::div(alpha*AIm)
//       - fvc::div(VReSource)
//     );
//     
//     
//     // Imaginary part
//     surfaceScalarField rAImAf
//     (
//         "rAImAf",
//         sigmaf*fvc::interpolate(1.0/AImEqn.A())
//     );
// 
//     surfaceScalarField VImSource
//     (
//         "VImSource",
//         rAImAf*(fvc::interpolate(fvc::grad(VIm, "grad(VImSource)")) & mesh.Sf())
//     );
//     
//     fvScalarMatrix VImEqn
//     (
//       - fvm::laplacian(rAImAf, VIm)
// //      ==
// //       + fvc::div(alpha*ARe)
//       - fvc::div(VImSource)
//     );

    // Real part
    fvScalarMatrix VReEqn
    (
      - fvm::laplacian(sigmaf, VRe)
//      ==
//       - fvc::div(alpha*AIm)
    );
    
    // Imaginary part
    fvScalarMatrix VImEqn
    (
      - fvm::laplacian(sigmaf, VIm)
//      ==
//       + fvc::div(alpha*ARe)
    );
