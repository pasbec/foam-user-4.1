
    // Source coefficient
    dimensionedScalar one
    (
        "one",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        1.0
    );

    // Real part
    fvScalarMatrix divAReEqn
    (
        fvm::Sp(one,divARe)
//       - fvc::div(ARe)
    );

    // Currently only work-around
    // since diagonal insertion
    // in block-matrix does not
    // work correctly
    divAReEqn.upper() = 0.0;

    // Imaginary part
    fvScalarMatrix divAImEqn
    (
        fvm::Sp(one,divAIm)
//       - fvc::div(AIm)
    );

    // Currently only work-around
    // since diagonal insertion
    // in block-matrix does not
    // work correctly
    divAImEqn.upper() = 0.0;
