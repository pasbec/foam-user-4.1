

    // Prepare block system
    fvBlockMatrix<vector10> AVEqn(AV);

    // Insert equations for A
    AVEqn.insertEquation(0, AReEqn);
    AVEqn.insertEquation(5, AImEqn);

    // Insert equations for G
    AVEqn.insertEquation(3, GReEqn);
    AVEqn.insertEquation(8, GImEqn);

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

//         // Coupling of real part of gauge potential equation
//         BlockLduSystem<vector, vector> GReInAIm(fvm::grad(GRe, "grad(GRe)")); GReInAIm *= alphai;
        BlockLduSystem<vector, scalar> AReInGRe(fvm::UDiv(ARe, "div(ARe)")); AReInGRe *= -1.0;
//         fvScalarMatrix mGReInVIm(fvm::laplacian(alphaf, GRe, "laplacian(alphaf,GRe)")); mGReInVIm *= -1.0;
//         BlockLduSystem<scalar, scalar> GReInVIm(mGReInVIm.mesh());
//         GReInVIm.diag() = mGReInVIm.diag();
//         GReInVIm.upper() = mGReInVIm.upper();

//         AVEqn.insertBlockCoupling(5, 3, GReInAIm, true);
        AVEqn.insertBlockCoupling(3, 0, AReInGRe, false);
//         AVEqn.insertBlockCoupling(9, 3, GReInVIm, true);

        // Coupling of imaginary part of gauge potential equation
//         BlockLduSystem<vector, vector> GImInARe(fvm::grad(GIm, "grad(GIm)")); GImInARe *= -alphai;
        BlockLduSystem<vector, scalar> AImInGIm(fvm::UDiv(AIm, "div(AIm)")); AImInGIm *= -1.0;
//         fvScalarMatrix mGImInVRe(fvm::laplacian(alphaf, GIm, "laplacian(alphaf,GIm)")); mGImInVRe *= 1.0;
//         BlockLduSystem<scalar, scalar> GImInVRe(mGImInVRe.mesh());
//         GImInVRe.diag() = mGImInVRe.diag();
//         GImInVRe.upper() = mGImInVRe.upper();
//
//         AVEqn.insertBlockCoupling(0, 8, GImInARe, true);
        AVEqn.insertBlockCoupling(8, 5, AImInGIm, false);
//         AVEqn.insertBlockCoupling(4, 8, GImInVRe, true);

        // Coupling of real part of electric scalar potential equation
        BlockLduSystem<vector, vector> VReInARe(fvm::grad(VRe, "grad(VRe)")); VReInARe *= sigmai;
        BlockLduSystem<vector, scalar> AImInVRe(fvm::UDiv(alphaf, AIm, "div(alphaf,AIm)")); AImInVRe *= 1.0;

        AVEqn.insertBlockCoupling(0, 4, VReInARe, true);
        AVEqn.insertBlockCoupling(4, 5, AImInVRe, false);

        // Coupling of imaginary part of electric scalar potential equation
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

    AVEqn.retrieveSolution(3, GRe.internalField());
    AVEqn.retrieveSolution(8, GIm.internalField());

    GRe.correctBoundaryConditions();
    GIm.correctBoundaryConditions();

    AVEqn.retrieveSolution(4, VRe.internalField());
    AVEqn.retrieveSolution(9, VIm.internalField());

    VRe.correctBoundaryConditions();
    VIm.correctBoundaryConditions();
