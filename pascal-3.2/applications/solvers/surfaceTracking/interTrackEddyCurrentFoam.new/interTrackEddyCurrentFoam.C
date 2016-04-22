/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "interTrackEddyCurrentApp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    using namespace Foam;

    argList::validOptions.insert("prefix", "name");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createRegionDynamicFvMeshUninitialized.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    interTrackEddyCurrentApp::Manager multiManager(args, runTime, regionMesh);

    interTrackApp::Manager& interTrackAppManager =
        multiManager.interTrackAppManager();

    eddyCurrentApp::Manager& eddyCurrentAppManager =
        multiManager.eddyCurrentAppManager();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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

// TODO: Magnetic update?

// TODO: Magnetic update?
        // Update mesh in buffer region
        if (true)
        {
            using namespace interTrackEddyCurrentApp;
            using namespace interTrackEddyCurrentApp::Region;

            // Calculate mesh velocity at fluid/buffer-interface
            // in buffer region from current boundary displacement
            // in fluid region
            regionMesh.patchMapMeshVelocityDirectMapped
            (
                Region::FLUID,
                Region::BUFFER
            );

            // Grab current points of buffer region as new point field
            pointField newPoints = regionMesh[Region::BUFFER].points();

            // Correct points for 2D-motion of buffer region
            twoDPointCorrector bTwoDPointCorr(regionMesh[Region::BUFFER]);
            bTwoDPointCorr.correctPoints(newPoints);

            // Use motionSolver to move and update mesh of buffer region
            regionMesh[Region::BUFFER].movePoints(newPoints);
            regionMesh[Region::BUFFER].update();
        }

// TODO: Magnetic update?
        // Update mesh of default region
        if (true)
        {
            using namespace interTrackEddyCurrentApp;
            using namespace interTrackEddyCurrentApp::Region;

            // Create new point field for default region
            // with current point positions of fluid region
            pointField newPoints = regionMesh.rmap(Region::FLUID);

            // Replace point positions of buffer region in
            // new point field for default region
            regionMesh.rmap(newPoints, Region::BUFFER);

            // Move and update mesh of default region
            regionMesh[Region::DEFAULT].movePoints(newPoints);
            regionMesh[Region::DEFAULT].update();
        }

// TODO: Magnetic update?
        // Update mesh of conductor region
        if (true)
        {
            using namespace interTrackEddyCurrentApp;
            using namespace interTrackEddyCurrentApp::Region;

            // Create new point field for conductor region
            // with current points of default region
            pointField newPoints = regionMesh.map(Region::CONDUCTOR);

            // Move and update mesh of conductor region
            regionMesh[Region::CONDUCTOR].movePoints(newPoints);
            regionMesh[Region::CONDUCTOR].update();
        }

// TODO: Magnetic update?
        // Solve eddy-current problem
        if (true)
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

// TODO: Magnetic update?
        if (true)
        {
            using namespace interTrackEddyCurrentApp;
            using namespace interTrackEddyCurrentApp::Region;

// TODO: Extrapolation
//             eddyCurrentAppManager.storage.FL().mapExtrapolate(Region::FLUID);
//             eddyCurrentAppManager.storage.pB().mapExtrapolate(Region::FLUID);
            eddyCurrentAppManager.storage().FL().mapCopyInternal(Region::FLUID);
            eddyCurrentAppManager.storage().pB().mapCopyInternal(Region::FLUID);
        }

// TODO: Magnetic update?
        // Set volume force in fluid region
        if (true)
        {
            using namespace interTrackApp;
            using namespace interTrackApp::Region;

            Manager& manager = interTrackAppManager;

            SM_GLOBALREGIONSCOPE(DEFAULT);

            volVectorField FL =
                eddyCurrentAppManager.storage().FL()[Region::DEFAULT];
            volScalarField pB =
                eddyCurrentAppManager.storage().pB()[Region::DEFAULT];

            // Use Lorentz force as volume force in fluid region
            storage.F() == FL;

// TODO: Magnetic pressure?
//             pB *= lorentzForceRotationalFactor;
//
//             F == FL + fvc::grad(pB);
//             F *= lorentzForceVolumeFactor;
        }

        // Solve fluid flow in fluid region
        {
            using namespace interTrackApp;
            using namespace interTrackApp::Region;

            Manager& manager = interTrackAppManager;

            SM_MANAGERSCOPE();
            SM_REGIONSCOPE(DEFAULT);

#           include "UpLoop.H"
        }

// TODO: TEST Buffer region?

// TODO: Write
    }

    return(0);
}


// ************************************************************************* //
