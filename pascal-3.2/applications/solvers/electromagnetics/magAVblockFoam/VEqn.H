
    // Real part
    surfaceScalarField rAReAf
    (
        "rAReAf",
        fvc::interpolate(1.0/AReEqn.A())
    );

    surfaceScalarField VReSource
    (
        "VReSource",
        rAReAf*(fvc::interpolate(sigma*fvc::grad(VRe, "grad(VReSource)")) & mesh.Sf())
    );
    
    fvScalarMatrix VReEqn
    (
      - fvm::laplacian(rAReAf, VRe)
//      ==
//         - fvc::div(alpha*AIm)
//         - fvc::div(VReSource)
    );
    
    
    // Imaginary part
    surfaceScalarField rAImAf
    (
        "rAImAf",
        fvc::interpolate(1.0/AImEqn.A())
    );

    surfaceScalarField VImSource
    (
        "VImSource",
        rAImAf*(fvc::interpolate(sigma*fvc::grad(VIm, "grad(VImSource)")) & mesh.Sf())
    );
    
    fvScalarMatrix VImEqn
    (
        - fvm::laplacian(rAImAf, VIm)
//      ==
//         - fvc::div(alpha*ARe)
//         - fvc::div(VImSource)
    );
