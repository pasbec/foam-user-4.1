/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    interTrackFoam

Description

\*---------------------------------------------------------------------------*/

#include "interTrackApp.H"

// TODO: Coupled solution of U and p?

// TODO: Read Zhang S., Zhao, X., General formulations for Rhie-Chow interpolation,
//       ASME Heat Transfer/Fluids Engineering Summer Conference, HT-FED04, Charlotte, USA, 2004

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    using namespace Foam;

    argList::validOptions.insert("prefix", "name");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    using namespace interTrackApp;
    using namespace interTrackApp::Region;

    Manager manager(args, runTime, mesh);

    manager.read();
    manager.init();

    SM_GLOBALREGIONSCOPE(DEFAULT);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    volScalarField& p = storage.p();
    volVectorField& U = storage.U();
    surfaceScalarField& phi = storage.phi();
    volScalarField& rho = storage.rho();
    volScalarField& mu = storage.mu();
    trackedSurface& interface = storage.interface();

// TODO: Two fluids
// #   include "setRefCell.H"

#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    while (manager.run())
    {
        interface.updateMesh(); // abs phi
        interface.updateDisplacementDirections();

        interface.predictPoints(); // abs phi

        // --- OUTER corrector loop
        while (control.loop())
        {
            interface.updateBoundaryConditions(); // rel phi

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phi,U);

            autoPtr<fvVector4Matrix> UpEqn;

#           include "UEqn.H"

            if (!settings.UpCoupled)
            {
                // --- Pressure corrector PISO loop
                while (control.correct())
                {
#                   include "pEqnSegregated.H"
                }
            }
            else
            {
#               include "pEqnCoupled.H"
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

#           include "freeSurfaceContinuityErrs.H"
        }

#       include "volContinuity.H"
    }

    return(0);
}

// ************************************************************************* //
