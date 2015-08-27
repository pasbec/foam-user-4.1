
//         scalarField one(mesh.size(), 1.0);
//         AVEqn.insertEquationCoupling(3, 3, -one);
//         AVEqn.insertEquationCoupling(8, 8, -one);

    dimensionedScalar one
    (
        "one",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        1.0
    );

    // Real part
    fvScalarMatrix divAReEqn
    (
      - fvm::Sp(one,divARe)
//       + fvc::div(ARe)
    );

    // Imaginary part
    fvScalarMatrix divAImEqn
    (
      - fvm::Sp(one,divAIm)
//       + fvc::div(AIm)
    );
