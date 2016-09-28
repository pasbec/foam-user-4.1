

    {
        if (Control::debug)
        {
            Info << "interTrackApp::Control | UpLoop.H : "
                << "Commencing INTERTRACK U-p loop."
                << endl;
        }

        // --- INTERTRACK corrector loop
        while (control.loop())
        {
            volScalarField& p = storage.p();
            volVectorField& U = storage.U();
            surfaceScalarField& phi = storage.phi();
            const volScalarField& rho = storage.rho();
            const volScalarField& mu = storage.mu();
            trackedSurface& interface = storage.interface();

            interface.updateBoundaryConditions(); // rel phi

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phi,U);

            autoPtr<fvVector4Matrix> UpEqn;

#           include "UpLoop_UEqn.H"

            if (!settings.UpCoupled)
            {
                // --- Pressure corrector PISO loop
                while (control.correct())
                {
#                   include "UpLoop_pEqnSegregated.H"
                }
            }
            else
            {
#               include "UpLoop_pEqnCoupled.H"
            }

            if
            (
                control.turbCorr()
            && storage.item_turbulence().enabled()
            )
            {
                storage.turbulence().correct();
            }

            // Make the fluxes absolute
            fvc::makeAbsolute(phi, U);

            interface.correctPoints(); // abs phi

// TODO: Continuity errors
// #           include "trackedSurfaceContinuityErrs.H"
        }
    }
