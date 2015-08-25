
    // ==================================================================//
    // Real part
    // ==================================================================//

    fvScalarMatrix VReEqn
    (
        fvm::laplacian(sigma, VRe)
     == fvc::div(alpha*AIm)
    );

    VReEqn.solve();

    // ==================================================================//
    // Imaginary part
    // ==================================================================//

    fvScalarMatrix VImEqn
    (
      - fvm::laplacian(sigma, VIm)
     == fvc::div(alpha*ARe)
    );

    VImEqn.solve();

    // ==================================================================//
    // Correct magnetic vector potential
    // ==================================================================//

    ARe += rOmega * fvc::grad(VIm);
    AIm -= rOmega * fvc::grad(VRe);
