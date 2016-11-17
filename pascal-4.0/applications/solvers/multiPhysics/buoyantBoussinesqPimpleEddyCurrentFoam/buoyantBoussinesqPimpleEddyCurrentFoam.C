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

// TODO: Wall functions for alphat?

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

    pimpleApp::Manager& pimpleAppManager =
        masterManager.pimpleAppManager();

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

        uniformDimensionedScalarField& omega0 = globalStorage.omega0();

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

        eddyCurrentAppManager.storage().F().mapExtrapolate(Region::FLUID);
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
            using namespace buoyantBoussinesqPimpleEddyCurrentApp;
            using namespace buoyantBoussinesqPimpleEddyCurrentApp::Region;

            Manager& manager = masterManager;

            SM_MANAGERSCOPE();

            if (Control::debug)
            {
                Info << "buoyantBoussinesqPimpleEddyCurrentApp::Control : "
                    << "Map/Extrapolate Lorentz-force in fluid region."
                    << endl;
            }

            eddyCurrentAppManager.storage().F().mapExtrapolate(Region::FLUID);
        }

        // Solve fluid flow
        {
            using namespace pimpleApp;
            using namespace pimpleApp::Region;

            Manager& manager = pimpleAppManager;

            SM_MANAGERSCOPE();
            SM_REGIONSCOPE(DEFAULT);

#           include "UpLoop.H"
        }
    }

    return(0);
}


// ************************************************************************* //
