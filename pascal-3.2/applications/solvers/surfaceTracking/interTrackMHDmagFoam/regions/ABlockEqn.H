
    // ==================================================================//
    // Omega and alpha
    // ==================================================================//

    dimensionedScalar omega
    (
        "omega",
        2.0 * mathematicalConstant::pi * frequency
    );

    volScalarField alpha = omega * sigma;

    // ==================================================================//
    // Block matrix preparations
    // ==================================================================//

    // Prepare block system
    BlockLduMatrix<vector6> AEqn(mesh);

    // Reference block diagonal and set it to zero
    Field<tensor6>& AEqnD = AEqn.diag().asSquare(); AEqnD = tensor6::zero;

    // Reference linear off-diagonal and set it to zero
    Field<vector6>& AEqnU = AEqn.upper().asLinear(); AEqnU = vector6::zero;
    Field<vector6>& AEqnL = AEqn.lower().asLinear(); AEqnL = vector6::zero;

    // Set up source
    vector6Field AEqnS (mesh.nCells(), vector6::zero);

    // ==================================================================//
    // Boundary conditions
    // ==================================================================//

    // Transfer the coupled interface list for processor/cyclic/etc. boundaries
    AEqn.interfaces() = A.boundaryField().blockInterfaces();

    // ==================================================================//
    // Real part
    // ==================================================================//

    fvVectorMatrix AReEqn
    (
        - physicalConstant::rMu0 * fvm::laplacian(ARe)
    //  - alpha * AIm // substituted by point-implicit block-coupling
        - jsRe
    );

    // ==================================================================//
    // Imaginary part
    // ==================================================================//

    fvVectorMatrix AImEqn
    (
        - physicalConstant::rMu0 * fvm::laplacian(AIm)
    //  + alpha * ARe // substituted by point-implicit block-coupling
        - jsIm
    );

    // ==================================================================//
    // Assemble segregated equations in block structure
    // ==================================================================//

    // Inset equations into block Matrix
    blockMatrixTools::insertEquation(0, AReEqn, AEqn, A.internalField(), AEqnS);
    blockMatrixTools::insertEquation(3, AImEqn, AEqn, A.internalField(), AEqnS);

    // ==================================================================//
    // Assemble point-implicit coupled equation part in block structure
    // ==================================================================//

    // Add off-diagonal terms
    forAll(AEqnD, i)
    {
        dimensionedScalar alphai = alpha[i];

        AEqnD[i](0,3) = -alphai.value() * mesh.V()[i];
        AEqnD[i](3,0) = +alphai.value() * mesh.V()[i];
        AEqnD[i](1,4) = -alphai.value() * mesh.V()[i];
        AEqnD[i](4,1) = +alphai.value() * mesh.V()[i];
        AEqnD[i](2,5) = -alphai.value() * mesh.V()[i];
        AEqnD[i](5,2) = +alphai.value() * mesh.V()[i];
    }

    // ==================================================================//
    // Solve the block matrix
    // ==================================================================//

    BlockSolverPerformance<vector6> solverPerf =
        BlockLduSolver<vector6>::New
        (
            A.name(),
            AEqn,
            mesh.solutionDict().solver(A.name())
        )->solve(A,AEqnS);

    solverPerf.print();

    // ==================================================================//
    // Retrieve the solution and update for next iteration
    // ==================================================================//

    blockMatrixTools::retrieveSolution(0, ARe.internalField(), A);
    blockMatrixTools::retrieveSolution(3, AIm.internalField(), A);

    ARe.correctBoundaryConditions();
    AIm.correctBoundaryConditions();
