
    {
        // Prepare block system
        fvBlockMatrix<vector8> AVEqn(AV);

        // Insert equations for A
        AVEqn.insertEquation(0, AReEqn);
        AVEqn.insertEquation(4, AImEqn);

        // Insert equations for V
        AVEqn.insertEquation(3, VReEqn);
        AVEqn.insertEquation(7, VImEqn);

        // References to interal fields
        const scalarField& sigmaIn = sigma.internalField();
        const scalarField& alphaIn = alpha.internalField();
        const vectorField& gradAlphaIn = gradAlpha.internalField();

        // Add coupling terms for ARe-AIm coupling
        AVEqn.insertEquationCoupling(0, 4,  scaleA*alphaIn);
        AVEqn.insertEquationCoupling(4, 0, -scaleA*alphaIn);
        AVEqn.insertEquationCoupling(1, 5,  scaleA*alphaIn);
        AVEqn.insertEquationCoupling(5, 1, -scaleA*alphaIn);
        AVEqn.insertEquationCoupling(2, 6,  scaleA*alphaIn);
        AVEqn.insertEquationCoupling(6, 2, -scaleA*alphaIn);

        // Coupling of real part of electric scalar potential equation
        BlockLduSystem<vector, vector> VReInARe(fvm::grad(VRe)); VReInARe *= -scaleA*sigmaIn;
//         BlockLduSystem<vector, scalar> AImInVRe(fvm::UDiv(sigmaf, AIm, "div(sigmaf,AIm)")); AImInVRe *= -scaleV*omega.value();

        AVEqn.insertBlockCoupling(0, 3, VReInARe, true);
//         AVEqn.insertBlockCoupling(3, 4, AImInVRe, false);
        AVEqn.insertEquationCoupling(3, 4, -scaleV*gradAlphaIn.component(0));
        AVEqn.insertEquationCoupling(3, 5, -scaleV*gradAlphaIn.component(1));
        AVEqn.insertEquationCoupling(3, 6, -scaleV*gradAlphaIn.component(2));

        // Coupling of imaginary part of electric scalar potential equation
        BlockLduSystem<vector, vector> VImInAIm(fvm::grad(VIm)); VImInAIm *= -scaleA*sigmaIn;
//         BlockLduSystem<vector, scalar> AReInVIm(fvm::UDiv(sigmaf, ARe, "div(sigmaf,ARe)")); AReInVIm *= scaleV*omega.value();

        AVEqn.insertBlockCoupling(4, 7, VImInAIm, true);
//         AVEqn.insertBlockCoupling(7, 0, AReInVIm, false);
        AVEqn.insertEquationCoupling(7, 0,  scaleV*gradAlphaIn.component(0));
        AVEqn.insertEquationCoupling(7, 1,  scaleV*gradAlphaIn.component(1));
        AVEqn.insertEquationCoupling(7, 2,  scaleV*gradAlphaIn.component(2));

        // Block coupled solver call
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
    }
