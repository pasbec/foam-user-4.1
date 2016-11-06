
    {
        volScalarField rAU("rAU", 1.0/UEqn.A());
        surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

        U = rAU*UEqn.H();

        surfaceScalarField phiU
        (
            (fvc::interpolate(U) & mesh.Sf()) + fvc::ddtPhiCorr(rAU, U, phi)
        );

        phi = phiU + rAUf*fvc::interpolate(rhok)*(g & mesh.Sf());

        while (control.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian(rAUf, p) == fvc::div(phi)
            );

// TODO: Read pRefCell and pRefValue with control
// #           include "setReference.H"
//             pEqn.setReference(pRefCell, pRefValue);

            pEqn.solve(mesh.solutionDict().solver(p.select(control.finalInnerIter())));

            if (control.finalNonOrthogonalIter())
            {
                phi -= pEqn.flux();
            }
        }

        U += rAU*fvc::reconstruct((phi - phiU)/rAUf);
        U.correctBoundaryConditions();

// TODO: Continuity errors
// #       include "continuityErrs.H"
    }
