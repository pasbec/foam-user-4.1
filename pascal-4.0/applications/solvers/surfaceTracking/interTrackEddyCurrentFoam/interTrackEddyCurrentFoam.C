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
    interTrackEddyCurrentFoam

Description

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#include "interTrackEddyCurrentApp.H"

// TODO: Fix meshPhi issues!

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    using namespace Foam;

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createRegionDynamicFvMeshUninitialized.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    interTrackEddyCurrentApp::Manager masterManager(args, runTime, regionMesh);

    eddyCurrentApp::Manager& eddyCurrentAppManager =
        masterManager.eddyCurrentAppManager();

    interTrackApp::Manager& interTrackAppManager =
        masterManager.interTrackAppManager();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    masterManager.read();
    masterManager.init();

// TODO: IDEAS? Buffer region may be out of sync with fluid region on write
// TODO: Make this nicer!
    // Init eddyCurrentApp
    {
        using namespace eddyCurrentApp;
        using namespace eddyCurrentApp::Region;

        Manager& manager = eddyCurrentAppManager;

        SM_GLOBALREGIONSCOPE(DEFAULT);

        uniformDimensionedScalarField& omega0 = globalStorage.omega0();

        {
#           include "AVLoop_AUpdate.H"

#           include "AVLoop_VGradUpdate.H"

#           include "BUpdate.H"

#           include "derivedFields.H"
        }
    }

// TODO: IDEAS? Buffer region may be out of sync with fluid region on write
// TODO: Make this nicer!
    // Init interTrackEddyCurrentApp
    {
        using namespace interTrackEddyCurrentApp;
        using namespace interTrackEddyCurrentApp::Region;

        Manager& manager = masterManager;

        SM_MANAGERSCOPE();

        eddyCurrentAppManager.storage().FL().mapExtrapolate(Region::FLUID);

// TODO: Magnetic pressure?
//         eddyCurrentAppManager.storage().pB().mapExtrapolate(Region::FLUID);
//
//         const volScalarField& pB =
//             eddyCurrentAppManager.storage().pB()[Region::FLUID];
//
//         pB *= lorentzForceRotationalFactor;
//
//         F = FL + fvc::grad(pB);
//         F *= lorentzForceVolumeFactor;
    }

    while (masterManager.run())
    {
        // Update mesh in fluid region and predict interface points
        {
            using namespace interTrackApp;
            using namespace interTrackApp::Region;

            Manager& manager = interTrackAppManager;

            SM_MANAGERSCOPE();
            SM_REGIONSCOPE(DEFAULT);

#           include "meshUpdate.H"
        }

        // Check for magnetic update
        Switch emUpdate = masterManager.control().loop();

        if (emUpdate)
        {
            using namespace interTrackEddyCurrentApp;
            using namespace interTrackEddyCurrentApp::Region;

            Manager& manager = masterManager;

            SM_MANAGERSCOPE();

            if (Control::debug)
            {
                Info << "interTrackEddyCurrentApp::Control : "
                    << "Update of electromagntic fields due."
                    << endl;
            }
        }

        // Update meshes in buffer, default and conducting region
        if (emUpdate)
        {
            using namespace interTrackEddyCurrentApp;
            using namespace interTrackEddyCurrentApp::Region;

            Manager& manager = masterManager;

            SM_MANAGERSCOPE();


            // Update mesh in buffer region
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if (Control::debug)
            {
                Info << "interTrackEddyCurrentApp::Control : "
                    << "Update mesh of buffer region."
                    << endl;
            }

            // Calculate mesh velocity at fluid/buffer-interface
            // in buffer region from current boundary displacement
            // difference
            manager.mesh().patchMapMeshVelocityDirectMapped
            (
                Region::FLUID,
                Region::BUFFER
            );

            // Use motionSolver to move and update mesh of buffer region
            manager.mesh()[Region::BUFFER].update();


            // Update mesh in default region
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if (Control::debug)
            {
                Info << "interTrackEddyCurrentApp::Control : "
                    << "Update mesh of default region."
                    << endl;
            }

            {
                // Create new point field for default region
                // with current point positions of fluid region
                pointField newPoints = manager.mesh().rmap(Region::FLUID);

                // Replace point positions of buffer region in
                // new point field for default region
                manager.mesh().rmap(newPoints, Region::BUFFER);

                // Move mesh points of default region
                manager.mesh()[Region::DEFAULT].movePoints(newPoints);
            }


            // Update mesh in conducting region
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if (Control::debug)
            {
                Info << "interTrackEddyCurrentApp::Control : "
                    << "Update mesh of conducting region."
                    << endl;
            }

            {
                // Create new point field for conductor region
                // with current points of default region
                pointField newPoints = manager.mesh().map(Region::CONDUCTOR);

                // Move mesh points of conductor region
                manager.mesh()[Region::CONDUCTOR].movePoints(newPoints);
            }
        }

        // Solve eddy-current problem
        if (emUpdate)
        {
            using namespace eddyCurrentApp;
            using namespace eddyCurrentApp::Region;

            Manager& manager = eddyCurrentAppManager;

            SM_GLOBALREGIONSCOPE(DEFAULT);

            uniformDimensionedScalarField& omega0 = globalStorage.omega0();

            {
#               include "A0BiotSavart.H"

#               include "AVLoop.H"

#               include "BUpdate.H"

#               include "derivedFields.H"
            }
        }

        // Map/Extrapolate and update volume force in fluid region
        if (emUpdate)
        {
            using namespace interTrackEddyCurrentApp;
            using namespace interTrackEddyCurrentApp::Region;

            Manager& manager = masterManager;

            SM_MANAGERSCOPE();

            if (Control::debug)
            {
                Info << "interTrackEddyCurrentApp::Control : "
                    << "Map/Extrapolate Lorentz-force in fluid region."
                    << endl;
            }

            eddyCurrentAppManager.storage().FL().mapExtrapolate(Region::FLUID);

// TODO: Magnetic pressure?
//             eddyCurrentAppManager.storage().pB().mapExtrapolate(Region::FLUID);
//
//             const volScalarField& pB =
//                 eddyCurrentAppManager.storage().pB()[Region::FLUID];
//
//             pB *= lorentzForceRotationalFactor;
//
//             F = FL + fvc::grad(pB);
//             F *= lorentzForceVolumeFactor;
        }

        // Solve fluid flow
        {
            using namespace interTrackApp;
            using namespace interTrackApp::Region;

            Manager& manager = interTrackAppManager;

            SM_MANAGERSCOPE();
            SM_REGIONSCOPE(DEFAULT);

#           include "UpLoop.H"
        }

// TODO: IDEAS? Buffer region may be out of sync with fluid region on write
// TODO: Write
    }

    return(0);
}


// ************************************************************************* //
