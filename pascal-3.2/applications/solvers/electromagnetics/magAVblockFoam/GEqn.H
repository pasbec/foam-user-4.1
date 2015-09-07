
    // Real part
    fvScalarMatrix GReEqn
    (
        Gscale * fvm::laplacian(GRe, "laplacian(GRe)")
//         Gscale * fvc::div(ARe)
    );



    // Imaginary part
    fvScalarMatrix GImEqn
    (
        Gscale * fvm::laplacian(GIm, "laplacian(GIm)")
//         Gscale * fvc::div(AIm)
    );
