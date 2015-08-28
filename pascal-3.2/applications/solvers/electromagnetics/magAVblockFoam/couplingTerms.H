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

    //
    BlockLduSystem<vector, vector> divAReInARe(fvm::grad(divARe)); divAReInARe *= rMu0.value()/divAfac;
    BlockLduSystem<vector, scalar> AReInDivARe(fvm::UDiv(ARe)); AReInDivARe *= -divAfac;

    //
    AVEqn.insertBlockCoupling(0, 3, divAReInARe, true);
    AVEqn.insertBlockCoupling(3, 0, AReInDivARe, false);

    //
    BlockLduSystem<vector, vector> divAImInAIm(fvm::grad(divAIm)); divAImInAIm *= rMu0.value()/divAfac;
    BlockLduSystem<vector, scalar> AImInDivAIm(fvm::UDiv(AIm)); AImInDivAIm *= -divAfac;

    //
    AVEqn.insertBlockCoupling(5, 8, divAImInAIm, true);
    AVEqn.insertBlockCoupling(8, 5, AImInDivAIm, false);

    //
    BlockLduSystem<vector, vector> VReInARe(fvm::grad(VRe)); VReInARe *= sigmai;
    BlockLduSystem<vector, scalar> AReInVIm(fvm::UDiv(alphaf,ARe));

    //
    AVEqn.insertBlockCoupling(0, 4, VReInARe, true);
    AVEqn.insertBlockCoupling(9, 0, AReInVIm, false);

    //
    BlockLduSystem<vector, vector> VImInAIm(fvm::grad(VIm)); VImInAIm *= sigmai;
    BlockLduSystem<vector, scalar> AImInVRe(fvm::UDiv(alphaf,AIm)); AImInVRe *= -1.0;

    //
    AVEqn.insertBlockCoupling(5, 9, VImInAIm, true);
    AVEqn.insertBlockCoupling(4, 5, AImInVRe, false);
}
