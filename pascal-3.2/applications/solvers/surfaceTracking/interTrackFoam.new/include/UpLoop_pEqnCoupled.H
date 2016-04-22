
    tmp<volScalarField> rAU(1.0/UEqn().A()); rAU().rename("rAU");
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));

// TODO: Boundary conditions for p?

    // Hint: phiP is a pressure source (Rhie-Chow?). It is important
    //       to note, that it already uses the boundary values of p!
    surfaceScalarField phiP
    (
        "phiP",
        rAUf*(fvc::interpolate(fvc::grad(p)) & mesh.Sf())
//         rAUf*(fvc::interpolate(fvc::grad(p, "grad(pSource)")) & mesh.Sf())
    );

    // The non-orthogonal corrections for the coupled approach
    // are completely realized within the outer loop
    {
        // Hint: The minus-signs in the pressure equation pEqn
        //       are necessary as this equation gets inserted
        //       into the block matrix UpEqn afterwards!
        fvScalarMatrix pEqn
        (
            -fvm::laplacian(rAUf, p) == -fvc::div(phiP)
        );

// TODO: Two fluids
// TODO: Read pRefCell and pRefValue with control
// #       include "setReference.H"
//         pEqn.setReference(pRefCell, pRefValue);

        UpEqn().insertEquation(3, pEqn);

        {
            BlockLduSystem<vector, vector> pInU(fvm::grad(p));
            BlockLduSystem<vector, scalar> UInp(fvm::UDiv(U));

            UpEqn().insertBlockCoupling(0, 3, pInU, true);
            UpEqn().insertBlockCoupling(3, 0, UInp, false);
        }

        UpEqn().solve();
        UpEqn().retrieveSolution(0, U.internalField());
        UpEqn().retrieveSolution(3, p.internalField());

        UpEqn.clear();

        U.correctBoundaryConditions();
        p.correctBoundaryConditions();

        // Calculate the new flux
        phi = (fvc::interpolate(U) & mesh.Sf()) + phiP + pEqn.flux();

        // Explicitly relax pressure for phiP
        p.relax();

        // Manually store residuls after first
        // non-orthogonal iteration
        control.storeResiduals(storage.Up().name());
    }

// TODO: Continuity errors
// #   include "continuityErrs.H"

// TODO: Port ddtCorr from OpenFOAM 3.0.x
//     {
//         Uf = fvc::interpolate(U);
//         surfaceVectorField n(mesh.Sf()/mesh.magSf());
//         Uf += n*(phi/mesh.magSf() - (n & Uf));
//     }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);
