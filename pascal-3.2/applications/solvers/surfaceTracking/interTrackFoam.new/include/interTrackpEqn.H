    {
        volScalarField rUA = 1.0/UEqn.A();

        U = rUA*UEqn.H();

        phi = (fvc::interpolate(U) & mesh.Sf());

// TODO: Port ddtCorr from OpenFOAM 3.0.x
//         bool ddtPhiCorr = false;
// #       include "calcPhi.H"
// #       include "correctPhiAtInterface.H"

// TODO: Two fluids
// #       include "scalePhi.H"

        // Non-orthogonal pressure corrector loop
        while (control.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian(rUA, p) == fvc::div(phi)
            );

#           include "setReference.H"
//             pEqn.setReference(pRefCell, pRefValue);

            pEqn.solve
            (
                mesh.solutionDict().solverDict
                (
                    p.select(control.finalInnerIter())
                )
            );

            if (control.finalNonOrthogonalIter())
            {
                phi -= pEqn.flux();
            }
        }

        p.relax();

#       include "continuityErrs.H"

        // Momentum corrector
        U -= rUA*fvc::grad(p);
        U.correctBoundaryConditions();
    }

