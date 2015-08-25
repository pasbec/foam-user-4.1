
    // Real part of magnetic vector potential
    fvVectorMatrix AReEqn
    (
        - rMu0 * fvm::laplacian(ARe)
    //  - alpha * AIm // substituted by point-implicit block-coupling
//         + sigma * fvc::grad(VRe)
        - jsRe
    );

    // Imaginary part of magnetic vector potential
    fvVectorMatrix AImEqn
    (
        - rMu0 * fvm::laplacian(AIm)
    //  + alpha * ARe // substituted by point-implicit block-coupling
//         + sigma * fvc::grad(VIm)
        - jsIm
    );
    
    // Prepare block system
    fvBlockMatrix<vector6> blockM(A);

    // Insert equations into block Matrix
    blockM.insertEquation(0, AReEqn);
    blockM.insertEquation(3, AImEqn);

    // Add off-diagonal coupling terms
    scalarField& coupling = alpha.internalField();

    blockM.insertEquationCoupling(0, 3, -coupling);
    blockM.insertEquationCoupling(3, 0,  coupling);
    blockM.insertEquationCoupling(1, 4, -coupling);
    blockM.insertEquationCoupling(4, 1,  coupling);
    blockM.insertEquationCoupling(2, 5, -coupling);
    blockM.insertEquationCoupling(5, 2,  coupling);

//     // Update source coupling: coupling terms eliminated from source
//     blockM.updateSourceCoupling();

    //- Block coupled solver call
    blockM.solve();

    // Retrieve solution
    blockM.retrieveSolution(0, ARe.internalField());
    blockM.retrieveSolution(3, AIm.internalField());

    ARe.correctBoundaryConditions();
    AIm.correctBoundaryConditions();
