

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
        const scalarField& sigmaIn = sigma.internalField();
        const scalarField& alphaIn = alpha.internalField();
        const vectorField& gradAlphaIn = gradAlpha.internalField();

        // Add coupling terms for ARe-AIm coupling
        AVEqn.insertEquationCoupling(0, 4,  Ascale*alphaIn);
        AVEqn.insertEquationCoupling(4, 0, -Ascale*alphaIn);
        AVEqn.insertEquationCoupling(1, 5,  Ascale*alphaIn);
        AVEqn.insertEquationCoupling(5, 1, -Ascale*alphaIn);
        AVEqn.insertEquationCoupling(2, 6,  Ascale*alphaIn);
        AVEqn.insertEquationCoupling(6, 2, -Ascale*alphaIn);

        // Coupling of real part of electric scalar potential equation
        BlockLduSystem<vector, vector> VReInARe(fvm::grad(VRe, "grad(VRe)")); VReInARe *= -Ascale*sigmaIn;
//         BlockLduSystem<vector, scalar> AImInVRe(fvm::UDiv(alphaf, AIm, "div(alphaf,AIm)")); AImInVRe *= -Vscale;

        AVEqn.insertBlockCoupling(0, 3, VReInARe, true);
//         AVEqn.insertBlockCoupling(3, 4, AImInVRe, false);
        AVEqn.insertEquationCoupling(3, 4, -Vscale*gradAlphaIn.component(0));
        AVEqn.insertEquationCoupling(3, 5, -Vscale*gradAlphaIn.component(1));
        AVEqn.insertEquationCoupling(3, 6, -Vscale*gradAlphaIn.component(2));

        // Coupling of imaginary part of electric scalar potential equation
        BlockLduSystem<vector, vector> VImInAIm(fvm::grad(VIm, "grad(VIm)")); VImInAIm *= -Ascale*sigmaIn;
//         BlockLduSystem<vector, scalar> AReInVIm(fvm::UDiv(alphaf, ARe, "div(alphaf,ARe)")); AReInVIm *= Vscale;

        AVEqn.insertBlockCoupling(4, 7, VImInAIm, true);
//         AVEqn.insertBlockCoupling(7, 0, AReInVIm, false);
        AVEqn.insertEquationCoupling(7, 0,  Vscale*gradAlphaIn.component(0));
        AVEqn.insertEquationCoupling(7, 1,  Vscale*gradAlphaIn.component(1));
        AVEqn.insertEquationCoupling(7, 2,  Vscale*gradAlphaIn.component(2));
    }

//     // Relax block-equation
//     Field<vector10>& AVx = AVEqn.psi().internalField();
//     Field<vector10>& AVb = AVEqn.source();
//     AVEqn.relax(AVx,AVb,0.7);

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
