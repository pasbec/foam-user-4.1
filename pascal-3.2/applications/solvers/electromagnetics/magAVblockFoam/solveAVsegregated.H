
    {
        // Add explicit parts of V (segregated coupling)
        fvVectorMatrix AReEqnPre = AReEqn - scaleA * sigma * fvc::grad(VRe);
        fvVectorMatrix AImEqnPre = AImEqn - scaleA * sigma * fvc::grad(VIm);

        // Prepare block system
        fvBlockMatrix<vector6> AEqn(A);

        // Insert equations for A
        AEqn.insertEquation(0, AReEqnPre);
        AEqn.insertEquation(3, AImEqnPre);

        // References to interal alpha field
        const scalarField& alphaIn = alpha.internalField();

        // Add coupling terms for ARe-AIm coupling
        AEqn.insertEquationCoupling(0, 3,  scaleA*alphaIn);
        AEqn.insertEquationCoupling(3, 0, -scaleA*alphaIn);
        AEqn.insertEquationCoupling(1, 4,  scaleA*alphaIn);
        AEqn.insertEquationCoupling(4, 1, -scaleA*alphaIn);
        AEqn.insertEquationCoupling(2, 5,  scaleA*alphaIn);
        AEqn.insertEquationCoupling(5, 2, -scaleA*alphaIn);

        // Block coupled solver call to solve A
        AEqn.solve();

        // Retrieve solution
        AEqn.retrieveSolution(0, ARe.internalField());
        AEqn.retrieveSolution(3, AIm.internalField());

        ARe.correctBoundaryConditions();
        AIm.correctBoundaryConditions();

        // Add explicit parts of A (segregated coupling)
        fvScalarMatrix VReEqnPre = VReEqn - scaleV * (gradAlpha&AIm);
        fvScalarMatrix VImEqnPre = VImEqn + scaleV * (gradAlpha&ARe);

        // Solve V
        VReEqnPre.solve();
        VImEqnPre.solve();

        VRe.correctBoundaryConditions();
        VIm.correctBoundaryConditions();
    }
