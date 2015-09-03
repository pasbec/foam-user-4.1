
    // Prepare block system
    fvBlockMatrix<vector8> AVEqn(AV);

    // Insert equations for A
    AVEqn.insertEquation(0, AReEqn);
    AVEqn.insertEquation(4, AImEqn);

    // Insert equations for V
    AVEqn.insertEquation(3, VReEqn);
    AVEqn.insertEquation(7, VImEqn);

    {
        // References to interal fields
        const scalarField& sigmai = sigma.internalField();
        const scalarField& alphai = alpha.internalField();

        // Add coupling terms for ARe-AIm coupling
        AVEqn.insertEquationCoupling(0, 4, -alphai);
        AVEqn.insertEquationCoupling(4, 0,  alphai);
        AVEqn.insertEquationCoupling(1, 5, -alphai);
        AVEqn.insertEquationCoupling(5, 1,  alphai);
        AVEqn.insertEquationCoupling(2, 6, -alphai);
        AVEqn.insertEquationCoupling(6, 2,  alphai);

        // Coupling of real part of scalar potential equation
        BlockLduSystem<vector, vector> VReInARe(fvm::grad(VRe, "grad(VRe)")); VReInARe *= sigmai;
        BlockLduSystem<vector, scalar> AImInVRe(fvm::UDiv(alphaf, AIm, "div(alphaf,AIm)")); AImInVRe *=1.0/Vscale;

        AVEqn.insertBlockCoupling(0, 3, VReInARe, true);
        AVEqn.insertBlockCoupling(3, 4, AImInVRe, false);

        // Coupling of imaginary part of scalar potential equation
        BlockLduSystem<vector, vector> VImInAIm(fvm::grad(VIm, "grad(VIm)")); VImInAIm *= sigmai;
        BlockLduSystem<vector, scalar> AReInVIm(fvm::UDiv(alphaf, ARe, "div(alphaf,ARe)")); AReInVIm *= -1.0/Vscale;

        AVEqn.insertBlockCoupling(4, 7, VImInAIm, true);
        AVEqn.insertBlockCoupling(7, 0, AReInVIm, false);
    }

    //- Block coupled solver call
    AVEqn.solve();

    // Retrieve solution
    AVEqn.retrieveSolution(0, ARe.internalField());
    AVEqn.retrieveSolution(4, AIm.internalField());

    ARe.correctBoundaryConditions();
    AIm.correctBoundaryConditions();

    AVEqn.retrieveSolution(3, VRe.internalField());
    AVEqn.retrieveSolution(7, VIm.internalField());

    VRe.correctBoundaryConditions();
    VIm.correctBoundaryConditions();
