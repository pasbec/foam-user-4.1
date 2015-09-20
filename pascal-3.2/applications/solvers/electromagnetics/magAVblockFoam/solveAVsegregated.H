
    {
        // Add explicit parts of V (segregated A-V coupling)
        fvVectorMatrix AReEqnSegregated = AReEqn - sigma * fvc::grad(VRe);
        fvVectorMatrix AImEqnSegregated = AImEqn - sigma * fvc::grad(VIm);

        // Prepare block system
        fvBlockMatrix<vector6> AEqn(A);

        // Lookup solver type for AEqn
        word solverType =
            mesh.solutionDict().subDict("solvers").
            subDict(AEqn.psi().name()).lookup("solver");

        // Force diag/upper/lower square for AMG (matrix for A is asymmetric)
        if (solverType == "AMG")
        {
            Field<tensor6>& AEqnD = AEqn.diag().asSquare(); AEqnD = tensor6::zero;
            Field<tensor6>& AEqnU = AEqn.upper().asSquare(); AEqnU = tensor6::zero;
            Field<tensor6>& AEqnL = AEqn.lower().asSquare(); AEqnL = tensor6::zero;
        }

        // Insert equations for A
        AEqn.insertEquation(0, AReEqnSegregated);
        AEqn.insertEquation(3, AImEqnSegregated);

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

//         // Add explicit parts of A (segregated A-V coupling)
//         fvScalarMatrix VReEqnSegregated = VReEqn;
//         fvScalarMatrix VImEqnSegregated = VImEqn;
//
//         bool useGradAlpha = true;
//
//         if (useGradAlpha)
//         {
//             VReEqnSegregated -= gradAlpha&AIm;
//             VImEqnSegregated += gradAlpha&ARe;
//         }
//         else
//         {
//             VReEqnSegregated -= fvc::div(alpha*AIm);
//             VImEqnSegregated += fvc::div(alpha*ARe);
//
//             VReEqnSegregated += alpha * fvc::div(AIm);
//             VImEqnSegregated -= alpha * fvc::div(ARe);
//         }
//
//         // Solve V
//         VReEqnSegregated.solve();
//         VImEqnSegregated.solve();
//
//         VRe.correctBoundaryConditions();
//         VIm.correctBoundaryConditions();
    }

    jRe ==   alpha * AIm - sigma * fvc::grad(VRe);
    jIm == - alpha * ARe - sigma * fvc::grad(VIm);

    BRe == fvc::curl(ARe);
    BIm == fvc::curl(AIm);
