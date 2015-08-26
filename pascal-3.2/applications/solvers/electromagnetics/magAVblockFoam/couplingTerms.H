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
    
    // 
    BlockLduSystem<vector, vector> VReInARe(fvm::grad(VRe)); VReInARe *= sigmai;
    BlockLduSystem<vector, scalar> AReInVIm(fvm::UDiv(sigmaf,ARe));

    // 
    AVEqn.insertBlockCoupling(0, 3, VReInARe, true);
    AVEqn.insertBlockCoupling(7, 0, AReInVIm, false);
    
    //
    BlockLduSystem<vector, vector> VImInAIm(fvm::grad(VIm)); VImInAIm *= sigmai;
    BlockLduSystem<vector, scalar> AImInVRe(fvm::UDiv(sigmaf,AIm));

    // 
    AVEqn.insertBlockCoupling(4, 7, VImInAIm, true);
    AVEqn.insertBlockCoupling(3, 4, AImInVRe, false);
}
