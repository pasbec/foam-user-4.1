
    // Source coefficient
    dimensionedScalar divAfacDim
    (
        "one",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        divAfac
    );

    // Real part
    fvScalarMatrix divAReEqn
    (
        fvm::Sp(divAfacDim,divARe)
//       - divFac * fvc::div(ARe)
    );

    // Currently only work-around
    // since diagonal insertion
    // in block-matrix does not
    // work correctly
    divAReEqn.upper() = 0.0;

    // Imaginary part
    fvScalarMatrix divAImEqn
    (
        fvm::Sp(divAfacDim,divAIm)
//       - divFac * fvc::div(AIm)
    );

    // Currently only work-around
    // since diagonal insertion
    // in block-matrix does not
    // work correctly
    divAImEqn.upper() = 0.0;
