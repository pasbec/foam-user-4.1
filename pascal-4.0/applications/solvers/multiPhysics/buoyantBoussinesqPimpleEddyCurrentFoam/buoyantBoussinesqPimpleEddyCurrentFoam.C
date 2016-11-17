/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    buoyantBoussinesqPimpleEddyCurrentFoam

Description

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#include "buoyantBoussinesqPimpleEddyCurrentApp.H"

// TODO: Pressure naming p vs. p_rgh (gh, ghf, hRef, ...)

// TODO: Relaxation of pressure?

// TODO: Pressure boundary conditions (F, bouyantPressure vs. fixedFluxPressure)

// TODO: Continuity errors

// TODO: Wall functions for kappat?

// TODO: Radiation (radiation.ST, rhoCpRef, ...)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    using namespace Foam;

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createRegionMeshUninitialized.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    buoyantBoussinesqPimpleEddyCurrentApp::Manager masterManager(args, runTime, regionMesh);

    eddyCurrentApp::Manager& eddyCurrentAppManager =
        masterManager.eddyCurrentAppManager();

    buoyantBoussinesqPimpleApp::Manager& buoyantBoussinesqPimpleAppManager =
        masterManager.buoyantBoussinesqPimpleAppManager();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    masterManager.read();
    masterManager.init();

// TODO: Make this nicer!
    // Init eddyCurrentApp
    {
        using namespace eddyCurrentApp;
        using namespace eddyCurrentApp::Region;

        Manager& manager = eddyCurrentAppManager;

        SM_GLOBALREGIONSCOPE(DEFAULT);

        uniformDimensionedScalarField& omega0 = storage.omega0();

        {
#           include "AVLoop_sigmaUpdate.H"
#           include "AVLoop_AUpdate.H"
#           include "AVLoop_VGradUpdate.H"

#           include "BUpdate.H"

#           include "derivedFields.H"
        }
    }

// TODO: Make this nicer!
    // Init buoyantBoussinesqPimpleEddyCurrentApp
    {
        using namespace buoyantBoussinesqPimpleEddyCurrentApp;
        using namespace buoyantBoussinesqPimpleEddyCurrentApp::Region;

        Manager& manager = masterManager;

        SM_MANAGERSCOPE();

        eddyCurrentAppManager.storage().F().rmap(Region::CONDUCTOR);
        eddyCurrentAppManager.storage().F().mapExtrapolate(Region::FLUID);

        eddyCurrentAppManager.storage().Q().rmap(Region::CONDUCTOR);
        eddyCurrentAppManager.storage().Q().mapExtrapolate(Region::THERMAL);
    }

    while (masterManager.run())
    {
        // Check for magnetic update
        Switch emUpdate =
            eddyCurrentAppManager.control().needsUpdate
            (
                buoyantBoussinesqPimpleEddyCurrentApp::Region::FLUID
            );

        if (emUpdate)
        {
            using namespace buoyantBoussinesqPimpleEddyCurrentApp;
            using namespace buoyantBoussinesqPimpleEddyCurrentApp::Region;

            Manager& manager = masterManager;

            SM_MANAGERSCOPE();

            if (Control::debug)
            {
                Info << "buoyantBoussinesqPimpleEddyCurrentApp::Control : "
                    << "Update of electromagntic fields due."
                    << endl;
            }
        }

        // Solve eddy-current problem
        if (emUpdate)
        {
            using namespace eddyCurrentApp;
            using namespace eddyCurrentApp::Region;

            Manager& manager = eddyCurrentAppManager;

            SM_GLOBALREGIONSCOPE(DEFAULT);

            uniformDimensionedScalarField& omega0 = storage.omega0();

            {
#               include "A0BiotSavart.H"

#               include "AVLoop.H"

#               include "BUpdate.H"

#               include "derivedFields.H"
            }
        }

        // Map/Extrapolate and update volume force in fluid region and
        // Joule heat in conducting region
        if (emUpdate)
        {
            using namespace buoyantBoussinesqPimpleEddyCurrentApp;
            using namespace buoyantBoussinesqPimpleEddyCurrentApp::Region;

            Manager& manager = masterManager;

            SM_MANAGERSCOPE();

            if (Control::debug)
            {
                Info << "buoyantBoussinesqPimpleEddyCurrentApp::Control : "
                    << "Map/Extrapolate Lorentz-force to fluid region."
                    << endl;
            }

            eddyCurrentAppManager.storage().F().rmap(Region::CONDUCTOR);
            eddyCurrentAppManager.storage().F().mapExtrapolate(Region::FLUID);

            if (Control::debug)
            {
                Info << "buoyantBoussinesqPimpleEddyCurrentApp::Control : "
                    << "Map/Extrapolate Joule-heat to thermal region."
                    << endl;
            }

            eddyCurrentAppManager.storage().Q().rmap(Region::CONDUCTOR);
            eddyCurrentAppManager.storage().Q().mapExtrapolate(Region::THERMAL);
        }

        // Solve thermal problem
        {
            using namespace buoyantBoussinesqPimpleEddyCurrentApp;
            using namespace buoyantBoussinesqPimpleEddyCurrentApp::Region;

            Manager& manager = masterManager;

            SM_MANAGERSCOPE();

            volScalarField& T = masterManager.storage().T()[Region::THERMAL];
            volScalarField& lambda = masterManager.storage().lambda()[Region::THERMAL];
            volScalarField& Q = eddyCurrentAppManager.storage().Q()[Region::THERMAL];

// TODO: rho
// TODO: Cp

// TODO: alphat & influence on lambda

// TODO
// //             alphat = turbulence.nut()/Prt;
// //             alphat.correctBoundaryConditions();
// //
// //             volScalarField alphaEff("alphaEff", turbulence.nu()/Pr + alphat);
//
//             fvScalarMatrix TEqn
//             (
//                 fvm::ddt(T)
//               + fvm::div(phi, T)
//               - fvm::laplacian(alphaEff, T)
// //              ==
// //                 radiation.ST(rhoCpRef, T)
//             );
//
//             TEqn.relax();
//
//             TEqn.solve();
//
// //             radiation.correct();
        }

        // Map/Extrapolate and update temperature in fluid region
        {
            using namespace buoyantBoussinesqPimpleEddyCurrentApp;
            using namespace buoyantBoussinesqPimpleEddyCurrentApp::Region;

            Manager& manager = masterManager;

            SM_MANAGERSCOPE();

            if (Control::debug)
            {
                Info << "buoyantBoussinesqPimpleEddyCurrentApp::Control : "
                    << "Map/Extrapolate temperature to fluid region."
                    << endl;
            }

            masterManager.storage().T().rmap(Region::THERMAL);
            masterManager.storage().T().mapExtrapolate(Region::THERMAL);
        }

        // Solve fluid flow
        {
            using namespace buoyantBoussinesqPimpleApp;
            using namespace buoyantBoussinesqPimpleApp::Region;

            Manager& manager = buoyantBoussinesqPimpleAppManager;

            SM_MANAGERSCOPE();
            SM_REGIONSCOPE(DEFAULT);

            if (Control::debug)
            {
                Info<< Control::typeName << " | UTpLoop.H : "
                    << "Commencing PIMPLE U-T-p loop."
                    << endl;
            }

            // --- PIMPLE corrector loop
            while (control.loop())
            {
                uniformDimensionedVectorField& g = storage.g();
                volScalarField& T = storage.T();
                volScalarField& p = storage.p();
                volVectorField& U = storage.U();
                surfaceScalarField& phi = storage.phi();
                uniformDimensionedScalarField& beta = storage.beta();
                uniformDimensionedScalarField& TRef = storage.TRef();
// TODO: Pr, Prt, alphat
//                 uniformDimensionedScalarField& Pr = storage.Pr();
//                 uniformDimensionedScalarField& Prt = storage.Prt();
                incompressible::turbulenceModel& turbulence = storage.turbulence();
                volScalarField& rhok = storage.rhok();

#               include "UTpLoop_UEqn.H"
#               include "UTpLoop_rhokUpdate.H"

                // --- Pressure corrector loop
                while (control.correct())
                {
#                   include "UTpLoop_pEqn.H"
                }

                if (control.turbCorr())
                {
                    storage.turbulence().correct();
                }

// TODO: Pr, Prt, alphat
// TODO: alphat & influence on lambda
            }
        }
    }

    return(0);
}


// ************************************************************************* //
