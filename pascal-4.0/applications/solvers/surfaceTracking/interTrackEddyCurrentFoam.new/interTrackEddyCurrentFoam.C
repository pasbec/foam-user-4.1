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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    using namespace Foam;

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createRegionDynamicFvMeshUninitialized.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    interTrackEddyCurrentApp::Manager multiManager(args, runTime, regionMesh);

    eddyCurrentApp::Manager& eddyCurrentAppManager =
        multiManager.eddyCurrentAppManager();

    interTrackApp::Manager& interTrackAppManager =
        multiManager.interTrackAppManager();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    multiManager.read();
    multiManager.init();

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

// TODO: Make this nicer!
    // Init interTrackEddyCurrentApp
    {
        using namespace interTrackEddyCurrentApp;
        using namespace interTrackEddyCurrentApp::Region;

        eddyCurrentAppManager.storage().FL().mapExtrapolate(Region::FLUID);

// TODO: Use pointer instead of copy!
        interTrackAppManager.regions().region_DEFAULT().storage().F() =
            eddyCurrentAppManager.storage().FL()[Region::FLUID];

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

    while (multiManager.run())
    {
        // Update mesh in fluid region and predict interface points
        {
            using namespace interTrackApp;
            using namespace interTrackApp::Region;

            Manager& manager = interTrackAppManager;

            SM_GLOBALREGIONSCOPE(DEFAULT);

#           include "meshUpdate.H"
        }

        // Check for magnetic update
        Switch emUpdate = multiManager.control().loop();

        // Update meshes in buffer, default and conducting region
        if (emUpdate)
        {
            using namespace interTrackEddyCurrentApp;
            using namespace interTrackEddyCurrentApp::Region;

            // Create new points
            pointField newPoints;

            // Update mesh in buffer region
            // ~~~~~~

            // Calculate mesh velocity at fluid/buffer-interface
            // in buffer region from current boundary displacement
            // in fluid region
            multiManager.mesh().patchMapMeshVelocityDirectMapped
            (
                Region::FLUID,
                Region::BUFFER
            );

            // Grab current points of buffer region as new point field
            newPoints = multiManager.mesh()[Region::BUFFER].points();

            // Correct points for 2D-motion of buffer region
            twoDPointCorrector bTwoDPointCorr(multiManager.mesh()[Region::BUFFER]);
            bTwoDPointCorr.correctPoints(newPoints);

            // Use motionSolver to move and update mesh of buffer region
            multiManager.mesh()[Region::BUFFER].movePoints(newPoints);
            multiManager.mesh()[Region::BUFFER].update();

            // Update mesh in default region
            // ~~~~~~

            // Create new point field for default region
            // with current point positions of fluid region
            newPoints = multiManager.mesh().rmap(Region::FLUID);

            // Replace point positions of buffer region in
            // new point field for default region
            multiManager.mesh().rmap(newPoints, Region::BUFFER);

            // Move and update mesh of default region
            multiManager.mesh()[Region::DEFAULT].movePoints(newPoints);
            multiManager.mesh()[Region::DEFAULT].update();

            // Update mesh in conducting region
            // ~~~~~~

            // Create new point field for conductor region
            // with current points of default region
            newPoints = multiManager.mesh().map(Region::CONDUCTOR);

            // Move and update mesh of conductor region
            multiManager.mesh()[Region::CONDUCTOR].movePoints(newPoints);
            multiManager.mesh()[Region::CONDUCTOR].update();
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

            eddyCurrentAppManager.storage().FL().mapExtrapolate(Region::FLUID);

// TODO: Use pointer instead of copy!
            interTrackAppManager.regions().region_DEFAULT().storage().F() =
                eddyCurrentAppManager.storage().FL()[Region::FLUID];

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

// TODO: Write
    }

    return(0);
}


// ************************************************************************* //
