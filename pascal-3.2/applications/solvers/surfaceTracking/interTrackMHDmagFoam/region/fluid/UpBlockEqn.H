
    volVector4Field& Up = *UpPtr;

    // ==================================================================//
    // Block matrix preparations
    // ==================================================================//

    // Create block matrix
    BlockLduMatrix<vector4> UpEqn(mesh);

    // Reference block diagonal and set it to zero
    Field<tensor4>& UpEqnD = UpEqn.diag().asSquare(); UpEqnD = tensor4::zero;

    // Reference linear off-diagonal and set it to zero
    Field<vector4>& UpEqnU = UpEqn.upper().asLinear(); UpEqnU = vector4::zero;
    Field<vector4>& UpEqnL = UpEqn.lower().asLinear(); UpEqnL = vector4::zero;

    // Set up source
    vector4Field UpEqnS (mesh.nCells(), vector4::zero);

    // ==================================================================//
    // Boundary conditions
    // ==================================================================//

    // Set block interfaces properly
    UpEqn.interfaces() = Up.boundaryField().blockInterfaces();

    // ==================================================================//
    // U Equation
    // ==================================================================//

    // Momentum equation
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U, "ddt(U)")
      + rho*fvm::div(phi, U, "div(phi,U)")
      + rho*turbulence->divDevReff(U)
      - F
    );

    if (outerCorr != nOuterCorr - 1)
    {
        Info << "Relaxing UEqn with factor: "
            << mesh.solutionDict().relaxationFactor("U") << nl;
        UEqn.relax(mesh.solutionDict().relaxationFactor("U"));
    }

    // ==================================================================//
    // p Equation
    // ==================================================================//

    volScalarField rAU ("rAU", 1.0/UEqn.A());

    // Pressure parts of the continuity equation
    surfaceScalarField rAUf
    (
        "rAUf",
        fvc::interpolate(rAU)
    );

    surfaceScalarField presSource
    (
        "presSource",
        rAUf * (fvc::interpolate(fvc::grad(p)) & mesh.Sf())
    );

    fvScalarMatrix pEqn
    (
      - fvm::laplacian(rAUf, p)
        ==
      - fvc::div(presSource)
    );

// FIXME [Low]: Parallel mode: Reference setting not yet implemented (only for twoFluids)
    // # include "setReference.H"

    // ==================================================================//
    // Assemble segregated equations in block structure
    // ==================================================================//

    blockMatrixTools::insertEquation(0, UEqn, UpEqn, Up.internalField(), UpEqnS);
    blockMatrixTools::insertEquation(3, pEqn, UpEqn, Up.internalField(), UpEqnS);

    // ==================================================================//
    // Assemble point-implicit coupled equation part in block structure
    // ==================================================================//

    // WARNING: It seems to be important that UInp and pInU reside inside
    //          a seperate scope to work correctly? Not sure, though!
    {
        blockVectorMatrix UInp(fvm::div(U,"blockDiv(U)"));
        blockVectorMatrix pInU(fvm::grad(p,"blockGrad(p)"));

        blockMatrixTools::insertBlockCoupling(3, 0, UInp, U, UpEqn, UpEqnS, false);
        blockMatrixTools::insertBlockCoupling(0, 3, pInU, p, UpEqn, UpEqnS, true);
    }

    // ==================================================================//
    // Solve the block matrix
    // ==================================================================//

    // Solve the block matrix
    BlockSolverPerformance<vector4> solverPerf =
        BlockLduSolver<vector4>::New
        (
            Up.name(),
            UpEqn,
            mesh.solutionDict().solver(Up.name())
        )->solve(Up, UpEqnS);

    solverPerf.print();

    // ==================================================================//
    // Retrieve the solution and update boundary conditions
    // ==================================================================//

    // Retrieve solution
    blockMatrixTools::retrieveSolution(0, U.internalField(), Up);
    blockMatrixTools::retrieveSolution(3, p.internalField(), Up);

    U.correctBoundaryConditions();
    p.correctBoundaryConditions();

    // ==================================================================//
    // Update phi, print continuity errors and relax pressure
    // ==================================================================//

    phi = (fvc::interpolate(U) & mesh.Sf()) + pEqn.flux() + presSource;

#   include "continuityErrs.H"

    // Explicitly relax p for momentum corrector
    if (outerCorr != nOuterCorr - 1)
    {
        scalar alpha = mesh.solutionDict().relaxationFactor("p");

        if (alpha < 1.0)
        {
            Info << "Relaxing p with factor: "
                << alpha << nl;

            p.relax(alpha);
        }
    }

    // ==================================================================//
    // Debugging stuff
    // ==================================================================//

    if (debug)
    {
        volScalarField& debug_divU = *debug_divUPtr;
        volScalarField& debug_divPhi = *debug_divPhiPtr;
        volVectorField& debug_Uphi = *debug_UphiPtr;
        volVectorField& debug_UmeshPhi = *debug_UmeshPhiPtr;

        scalar meshH = pow(max(mesh.V()),1.0/3.0).value();

        debug_divU = fvc::div(U) * meshH;
        debug_divPhi = fvc::div(phi) * meshH;

        debug_Uphi = fvc::reconstruct(phi);
        debug_UmeshPhi = fvc::reconstruct(fvc::meshPhi(U));
    }