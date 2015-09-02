

    // Prepare block system
    fvBlockMatrix<vector10> AVEqn(AV);

    // Insert equations for A
    AVEqn.insertEquation(0, AReEqn);
    AVEqn.insertEquation(5, AImEqn);

    // Insert equations for divA
    AVEqn.insertEquation(3, divAReEqn);
    AVEqn.insertEquation(8, divAImEqn);

    // Insert equations for V
    AVEqn.insertEquation(4, VReEqn);
    AVEqn.insertEquation(9, VImEqn);

    {
        // References to interal fields
        const scalarField& sigmai = sigma.internalField();
        const scalarField& alphai = alpha.internalField();

        // Add coupling terms for ARe-AIm coupling
        AVEqn.insertEquationCoupling(0, 5, -alphai);
        AVEqn.insertEquationCoupling(5, 0,  alphai);
        AVEqn.insertEquationCoupling(1, 6, -alphai);
        AVEqn.insertEquationCoupling(6, 1,  alphai);
        AVEqn.insertEquationCoupling(2, 7, -alphai);
        AVEqn.insertEquationCoupling(7, 2,  alphai);

        // Gauge-coupling for real part of A
        BlockLduSystem<vector, vector> divAReInARe(fvm::grad(divARe, "grad(divARe)")); divAReInARe *= rMu0.value();
        BlockLduSystem<vector, scalar> AReInDivARe(fvm::UDiv(ARe, "div(ARe)")); AReInDivARe *= -1.0;

        AVEqn.insertBlockCoupling(0, 3, divAReInARe, true);
        AVEqn.insertBlockCoupling(3, 0, AReInDivARe, false);

        // Gauge-coupling for imaginary part of A
        BlockLduSystem<vector, vector> divAImInAIm(fvm::grad(divAIm, "grad(divAIm)")); divAImInAIm *= rMu0.value();
        BlockLduSystem<vector, scalar> AImInDivAIm(fvm::UDiv(AIm, "div(AIm)")); AImInDivAIm *= -1.0;

        AVEqn.insertBlockCoupling(5, 8, divAImInAIm, true);
        AVEqn.insertBlockCoupling(8, 5, AImInDivAIm, false);

        // Coupling of real part of scalar potential equation
        BlockLduSystem<vector, vector> VReInARe(fvm::grad(VRe, "grad(VRe)")); VReInARe *= sigmai;
        BlockLduSystem<vector, scalar> AImInVRe(fvm::UDiv(alphaf, AIm, "div(alphaf,AIm)")); AImInVRe *= 1.0;

        AVEqn.insertBlockCoupling(0, 4, VReInARe, true);
        AVEqn.insertBlockCoupling(4, 5, AImInVRe, false);

        // Coupling of imaginary part of scalar potential equation
        BlockLduSystem<vector, vector> VImInAIm(fvm::grad(VIm, "grad(VIm)")); VImInAIm *= sigmai;
        BlockLduSystem<vector, scalar> AReInVIm(fvm::UDiv(alphaf, ARe, "div(alphaf,ARe)")); AReInVIm *= -1.0;

        AVEqn.insertBlockCoupling(5, 9, VImInAIm, true);
        AVEqn.insertBlockCoupling(9, 0, AReInVIm, false);
    }

    //- Block coupled solver call
    AVEqn.solve();

    // Retrieve solution
    AVEqn.retrieveSolution(0, ARe.internalField());
    AVEqn.retrieveSolution(5, AIm.internalField());

    ARe.correctBoundaryConditions();
    AIm.correctBoundaryConditions();

    AVEqn.retrieveSolution(3, divARe.internalField());
    AVEqn.retrieveSolution(8, divAIm.internalField());

    divARe.correctBoundaryConditions();
    divAIm.correctBoundaryConditions();

    AVEqn.retrieveSolution(4, VRe.internalField());
    AVEqn.retrieveSolution(9, VIm.internalField());

    VRe.correctBoundaryConditions();
    VIm.correctBoundaryConditions();
