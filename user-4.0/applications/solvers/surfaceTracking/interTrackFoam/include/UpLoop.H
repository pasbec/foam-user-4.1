

    {
        if (Control::debug)
        {
            Info<< Control::typeName << " | UpLoop.H : "
                << "Commencing PIMPLE U-p loop."
                << endl;
        }

        // --- PIMPLE corrector loop
        while (control.loop())
        {
            uniformDimensionedVectorField& g = storage.g();
            volScalarField& p = storage.p();
            volVectorField& U = storage.U();
            surfaceScalarField& phi = storage.phi();
// TODO: Make trackedSurface use real rho and mu fields!
            volScalarField& rho = storage.rho();
            twoPhaseMixture& transport = storage.transport();
            incompressible::turbulenceModel& turbulence = storage.turbulence();
            trackedSurface& interface = storage.interface();

            if (settings.relToUinf)
            {
// TODO: Real forcing!
                vector Uinf = storage.Uinf().value();

                storage.finf().value() =
                    (interface.totalNormalForce()
                  + interface.totalTangentialForce())
                  & Uinf/(mag(Uinf) + VSMALL);

                scalar dUinf =
                    storage.dUinfRelax().value()
                  * storage.finf().value()
                  / max(rho*transport.nu()).value()
                  / runTime.deltaT().value();

                // Update free-stream velocity Uinf
                storage.Uinf().value() = Uinf + dUinf*Uinf/(mag(Uinf) + VSMALL);
// TODO: Real forcing!

                Info<< "Free-stream velocity adjustment:"
                    << " finf = " << storage.finf().value()
                    << " , Uinf = " << storage.Uinf().value()
                    << " , dUinfRelax = " << storage.dUinfRelax().value()
                    << " , deltaT = " << runTime.deltaT().value()
                    << " , dUinf = " << dUinf
                    << endl;

                // Remove old counterflow velocity U0 from U
                {
                    // WARNING: FORCE assign boundary values!
                    U == U - storage.U0();
                    U.correctBoundaryConditions();

                    phi -= storage.phi0();
                }

                // Scale counterflow velocity according to new free-stream velocity
                // and sdjust boundary conditions
                {
                    storage.U0() *= mag(storage.Uinf().value())/(mag(Uinf) + VSMALL);
                    storage.phi0() *= mag(storage.Uinf().value())/(mag(Uinf) + VSMALL);

                    forAll(storage.U0().boundaryField(), patchI)
                    {
                        if
                        (
                            storage.U0().boundaryField()[patchI].type()
                        == fixedValueFvPatchVectorField::typeName
                        && patchI != interface.aPatchID()
                        )
                        {
                            storage.U0().boundaryField()[patchI] == storage.Uinf().value();
                        }
                    }

                    storage.U0().correctBoundaryConditions();
                }

                // Update counterflow velocity U0
                int maxIter = 10;
                scalar tolerance = 1e-3;
                scalar residual = VGREAT;
                for (int iter=0; iter<maxIter; iter++)
                {
                    volScalarField& p0 = storage.p0();
                    volVectorField& U0 = storage.U0();
                    surfaceScalarField& phi0 = storage.phi0();

                    volScalarField mu
                    (
                        "mu",
                        rho*transport.nu()
                    );

                    surfaceScalarField muf
                    (
                        "muf",
                        fvc::interpolate(mu)
                    );

                    surfaceScalarField rhoPhi0
                    (
                        "rho*phi0",
                        fvc::interpolate(rho)*phi0
                    );

                    fvVectorMatrix U0Eqn
                    (
                        fvm::div(rhoPhi0, U0)
                      - fvm::laplacian(muf, U0)
                      - (fvc::grad(U0) & fvc::grad(muf))
                    );

                    U0Eqn.relax();

                    solve(U0Eqn == -fvc::grad(p0));

                    {
                        volScalarField rAU("rAU", 1.0/U0Eqn.A());
                        surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

                        U0 = rAU*U0Eqn.H();

                        phi0 = fvc::interpolate(U0) & mesh.Sf();

                        // Non-orthogonal pressure corrector loop
                        while (control.correctNonOrthogonal())
                        {
                            fvScalarMatrix p0Eqn
                            (
                                fvm::laplacian(rAUf,p0) == fvc::div(phi0)
                            );

                            p0Eqn.setReference(control.pRefCell(), control.pRefValue());

                            residual = p0Eqn.solve().initialResidual();

                            if (control.finalNonOrthogonalIter())
                            {
                                phi0 -= p0Eqn.flux();
                            }
                        }

                        p0.relax();

                        U0 -= rAU*fvc::grad(p0);
                        U0.correctBoundaryConditions();
                    }

                    if (residual < tolerance)
                    {
                        break;
                    }
                }

                // Add new counterflow velocity U0 to U
                {
                    // WARNING: FORCE assign boundary values!
                    U == U + storage.U0();
                    U.correctBoundaryConditions();

                    phi += storage.phi0();
                }
            }

            interface.updateBoundaryConditions(); // abs phi

            surfaceScalarField gxf("gxf", g & mesh.Cf());

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phi, U);

#           include "UpLoop_UEqn.H"

            // --- Pressure corrector PISO loop
            while (control.correct())
            {
#               include "UpLoop_pEqn.H"
            }

            if (control.turbCorr())
            {
                turbulence.correct();
            }

#           include "UpLoop_cEqn.H"
#           include "UpLoop_rhoUpdate.H"

            // Make the fluxes absolute
            fvc::makeAbsolute(phi, U);

            interface.updatePoints(); // abs phi

// TODO: Continuity errors
// #           include "trackedSurfaceContinuityErrs.H"
        }
    }
