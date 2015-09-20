
    {
        // Prepare block system
        fvBlockMatrix<vector8> AVEqn(AV);

        // Insert equations for A
        AVEqn.insertEquation(0, AReEqn);
        AVEqn.insertEquation(4, AImEqn);

        // Insert equations for V
        AVEqn.insertEquation(3, VReEqn);
        AVEqn.insertEquation(7, VImEqn);

        // References to internal fields
        const scalarField& sigmaIn = sigma.internalField();
        const scalarField& alphaIn = alpha.internalField();
        const vectorField& gradAlphaIn = gradAlpha.internalField();

        // Add coupling terms for ARe-AIm coupling
        AVEqn.insertEquationCoupling(0, 4,  alphaIn);
        AVEqn.insertEquationCoupling(4, 0, -alphaIn);
        AVEqn.insertEquationCoupling(1, 5,  alphaIn);
        AVEqn.insertEquationCoupling(5, 1, -alphaIn);
        AVEqn.insertEquationCoupling(2, 6,  alphaIn);
        AVEqn.insertEquationCoupling(6, 2, -alphaIn);

        // Coupling part of electric scalar potential equation in A
        BlockLduSystem<vector, vector> VReInARe(fvm::grad(VRe)); VReInARe *= -sigmaIn;
        BlockLduSystem<vector, vector> VImInAIm(fvm::grad(VIm)); VImInAIm *= -sigmaIn;

        AVEqn.insertBlockCoupling(0, 3, VReInARe, true);
        AVEqn.insertBlockCoupling(4, 7, VImInAIm, true);

        // Coupling part of electric scalar potential equation from A
        bool useGradAlpha = true;

        if (useGradAlpha)
        {
            AVEqn.insertEquationCoupling(3, 4, -gradAlphaIn.component(0));
            AVEqn.insertEquationCoupling(3, 5, -gradAlphaIn.component(1));
            AVEqn.insertEquationCoupling(3, 6, -gradAlphaIn.component(2));

            AVEqn.insertEquationCoupling(7, 0,  gradAlphaIn.component(0));
            AVEqn.insertEquationCoupling(7, 1,  gradAlphaIn.component(1));
            AVEqn.insertEquationCoupling(7, 2,  gradAlphaIn.component(2));
        }
        else
        {
            BlockLduSystem<vector, scalar> AImInVRe(fvm::UDiv(sigmaf, AIm, "div(sigmaf,AIm)"));
            BlockLduSystem<vector, scalar> AReInVIm(fvm::UDiv(sigmaf, ARe, "div(sigmaf,ARe)"));

            AImInVRe *= -omega.value();
            AReInVIm *=  omega.value();

            BlockLduSystem<vector, scalar> AImInVReDiv(fvm::UDiv(AIm));
            BlockLduSystem<vector, scalar> AReInVImDiv(fvm::UDiv(ARe));

            AImInVReDiv *= alphaIn;
            AReInVImDiv *= alphaIn;

            AImInVRe += AImInVReDiv;
            AReInVIm -= AReInVImDiv;

            AVEqn.insertBlockCoupling(3, 4, AImInVRe, false);
            AVEqn.insertBlockCoupling(7, 0, AReInVIm, false);
        }

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

    jRe ==   alpha * AIm - sigma * fvc::grad(VRe);
    jIm == - alpha * ARe - sigma * fvc::grad(VIm);

    BRe == fvc::curl(ARe);
    BIm == fvc::curl(AIm);
