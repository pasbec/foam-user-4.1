
    rAU = 1.0/UEqn.A();

    // Pressure parts of the continuity equation
    surfaceScalarField rAUf
    (
        "rAUf",
        fvc::interpolate(rAU)
    );

    surfaceScalarField presSource
    (
        "presSource",
        rAUf*(fvc::interpolate(fvc::grad(p)) & mesh.Sf())
    );

    fvScalarMatrix pEqn
    (
    - fvm::laplacian(rAUf, p)
    ==
    - fvc::div(presSource)
    );

    // FIXME [Low]: Parallel mode: Reference setting not yet implemented (only for twoFluids)
    // # include "setReference.H"

    if (outerCorr != nOuterCorr - 1)
    {
        Info << "Relaxing pEqn with factor: " << mesh.solutionDict().relaxationFactor("p") << nl;
        pEqn.relax(mesh.solutionDict().relaxationFactor("p"));
    }

    blockVectorMatrix UInp(fvm::div(U,"BLOCK[fvm::div(U)]"));