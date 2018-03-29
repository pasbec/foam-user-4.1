

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
                volScalarField mu
                (
                    "mu",
                    rho*transport.nu()
                );

                vector oldUinf = storage.Uinf().value();

                scalar finf =
                     (interface.totalNormalForce()
                   + interface.totalTangentialForce())
                   & oldUinf/mag(oldUinf);

                scalar dUinf =
                    storage.dUinfRelax().value() * finf
                  / (max(mu).value() * mag(oldUinf) * runTime.deltaT().value());

                vector Uinf = oldUinf + dUinf*oldUinf/mag(oldUinf);

                Info<< "Free stream adjustment:"
                    << " finf = " << finf
                    << " , Uinf = " << Uinf
                    << endl;

                storage.finf().value() = finf;
                storage.Uinf().value() = Uinf;

                forAll(U.boundaryField(), patchI)
                {
                    if
                    (
                        U.boundaryField()[patchI].type()
                        == fixedValueFvPatchVectorField::typeName
                    )
                    {
                        U.boundaryField()[patchI] == Uinf;
                    }
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

            if (settings.relToUinf)
            {
                phi -= storage.Uinf() & mesh.Sf();
            }

#           include "UpLoop_cEqn.H"
#           include "UpLoop_rhoUpdate.H"

            if (settings.relToUinf)
            {
                phi += storage.Uinf() & mesh.Sf();
            }

            // Make the fluxes absolute
            fvc::makeAbsolute(phi, U);

            interface.updatePoints(); // abs phi

// TODO: Continuity errors
// #           include "trackedSurfaceContinuityErrs.H"
        }
    }
