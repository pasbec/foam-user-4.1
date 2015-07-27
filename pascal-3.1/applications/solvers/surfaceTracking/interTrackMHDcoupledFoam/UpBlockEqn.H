
    // Insert U and p equations
    blockMatrixTools::insertEquation(0, UEqn, A, x, b);
    blockMatrixTools::insertEquation(3, pEqn, A, x, b);

    // Insert U-p coupling
    blockMatrixTools::insertBlockCoupling(3, 0, UInp, U, A, b, false);
    blockMatrixTools::insertBlockCoupling(0, 3, pInU, p, A, b, true);

    // Solve the block matrix
    BlockSolverPerformance<vector4> solverPerf =
        BlockLduSolver<vector4>::New
        (
            word("Up"),
            A,
            mesh.solutionDict().solver("Up")
        )->solve(Up, b);

    solverPerf.print();

    // Retrieve solution
    blockMatrixTools::retrieveSolution(0, U.internalField(), Up);
    blockMatrixTools::retrieveSolution(3, p.internalField(), Up);

    // Remove swirl
    if (removeSwirl > -1)
    {
        U.field().replace(removeSwirl, 0.0);
        U.boundaryField().replace(removeSwirl, 0.0);
    }

    U.correctBoundaryConditions();
    p.correctBoundaryConditions();

    phi = (fvc::interpolate(U) & mesh.Sf()) + pEqn.flux() + presSource;

#   include "continuityErrs.H"