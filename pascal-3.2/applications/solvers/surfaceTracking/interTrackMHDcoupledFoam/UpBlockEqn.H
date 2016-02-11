
    // ==================================================================//
    // Block matrix preparations
    // ==================================================================//

    // Initialize the Up block system (matrix, source and reference to Up)
    fvBlockMatrix<vector4> UpEqn(Up);

    // ==================================================================//
    // U Equation
    // ==================================================================//

    // Momentum equation
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U, "ddt(U)")
      + rho*fvm::div(phi, U, "div(phi,U)")
      + rho*turbulence->divDevReff(U)
      - FL
    );

    UEqn.relax();

    if (outerCorr != nOuterCorr - 1)
    {
        Info << "Relaxing UEqn with factor: "
            << mesh.solutionDict().relaxationFactor("U") << nl;
        UEqn.relax(mesh.solutionDict().relaxationFactor("U"));
    }

    UpEqn.insertEquation(0, UEqn);

    // ==================================================================//
    // p Equation
    // ==================================================================//

    // Pressure parts of the continuity equation
    surfaceScalarField rUAf
    (
        "rUAf",
        fvc::interpolate(1.0/UEqn.A())
    );

    surfaceScalarField presSource
    (
        "presSource",
        rUAf*(fvc::interpolate(fvc::grad(p, "grad(pSource)")) & mesh.Sf())
    );

    fvScalarMatrix pEqn
    (
      - fvm::laplacian(rUAf, p)
     ==
      - fvc::div(presSource)
    );

// FIXME [Low]: Parallel mode: Reference setting not yet implemented (only for twoFluids)
    // # include "setReference.H"

    UpEqn.insertEquation(3, pEqn);

    // ==================================================================//
    // Assemble coupling
    // ==================================================================//

    // WARNING: It seems to be important that UInp and pInU reside inside
    //          a seperate scope to work correctly? Not sure, though!
    {
        // Calculate grad p coupling matrix. Needs to be here if one uses
        // gradient schemes with limiters.  VV, 9/June/2014
        BlockLduSystem<vector, vector> pInU(fvm::grad(p));

        // Calculate div U coupling.  Could be calculated only once since
        // it is only geometry dependent.  VV, 9/June/2014
        BlockLduSystem<vector, scalar> UInp(fvm::UDiv(U));

        // Last argument in insertBlockCoupling says if the column direction
        // should be incremented. This is needed for arbitrary positioning
        // of U and p in the system. This could be better. VV, 30/April/2014
        UpEqn.insertBlockCoupling(0, 3, pInU, true);
        UpEqn.insertBlockCoupling(3, 0, UInp, false);
    }

    // ==================================================================//
    // Solve the block matrix
    // ==================================================================//

    UpEqn.solve();

    // ==================================================================//
    // Retrieve the solution and update boundary conditions
    // ==================================================================//

    // Retrieve solution
    UpEqn.retrieveSolution(0, U.internalField());
    UpEqn.retrieveSolution(3, p.internalField());

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
