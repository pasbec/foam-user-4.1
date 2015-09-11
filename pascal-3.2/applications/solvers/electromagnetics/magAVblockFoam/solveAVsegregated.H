
    {
        // Add explicit parts of V (segregated coupling)
        fvVectorMatrix AReEqnPre = AReEqn - sigma * fvc::grad(VRe);
        fvVectorMatrix AImEqnPre = AImEqn - sigma * fvc::grad(VIm);

        // Prepare block system
        fvBlockMatrix<vector6> AEqn(A);

        // Lookup solver type for AEqn
        word solverType =
            mesh.solutionDict().subDict("solvers").
            subDict(AEqn.psi().name()).lookup("solver");

        // Force diag/upper/lower square for AMG (matrix is asymmetric)
        if (solverType == "AMG")
        {
            Field<tensor6>& AEqnD = AEqn.diag().asSquare(); AEqnD = tensor6::zero;
            Field<tensor6>& AEqnU = AEqn.upper().asSquare(); AEqnU = tensor6::zero;
            Field<tensor6>& AEqnL = AEqn.lower().asSquare(); AEqnL = tensor6::zero;
        }

        // Insert equations for A
        AEqn.insertEquation(0, AReEqnPre);
        AEqn.insertEquation(3, AImEqnPre);

        // References to internal alpha field
        const scalarField& alphaIn = alpha.internalField();

        // Add coupling terms for ARe-AIm coupling
        AEqn.insertEquationCoupling(0, 3,  alphaIn);
        AEqn.insertEquationCoupling(3, 0, -alphaIn);
        AEqn.insertEquationCoupling(1, 4,  alphaIn);
        AEqn.insertEquationCoupling(4, 1, -alphaIn);
        AEqn.insertEquationCoupling(2, 5,  alphaIn);
        AEqn.insertEquationCoupling(5, 2, -alphaIn);

        // Block coupled solver call to solve A
        AEqn.solve();

        // Retrieve solution
        AEqn.retrieveSolution(0, ARe.internalField());
        AEqn.retrieveSolution(3, AIm.internalField());

        ARe.correctBoundaryConditions();
        AIm.correctBoundaryConditions();

        // Add explicit parts of A (segregated coupling)
        fvScalarMatrix VReEqnPre = VReEqn - (gradAlpha&AIm);
        fvScalarMatrix VImEqnPre = VImEqn + (gradAlpha&ARe);

        // Solve V
        VReEqnPre.solve();
        VImEqnPre.solve();

        VRe.correctBoundaryConditions();
        VIm.correctBoundaryConditions();
    }
