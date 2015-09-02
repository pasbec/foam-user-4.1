
    // Real part
    fvScalarMatrix divAReEqn
    (
        fvm::Sp(scalar(1.0),divARe)
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
        fvm::Sp(scalar(1.0),divAIm)
//       - fvc::div(AIm)
    );

    // Currently only work-around
    // since diagonal insertion
    // in block-matrix does not
    // work correctly
    divAImEqn.upper() = 0.0;
