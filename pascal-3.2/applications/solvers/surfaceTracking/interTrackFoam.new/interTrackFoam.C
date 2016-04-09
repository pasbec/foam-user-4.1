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

#include "interTrackManager.H"

#ifndef namespaceFoam
#define namespaceFoam
    using namespace Foam;
#endif

// TODO: Coupled solution of U and p?

// TODO: Read Zhang S., Zhao, X., General formulations for Rhie-Chow interpolation,
//       ASME Heat Transfer/Fluids Engineering Summer Conference, HT-FED04, Charlotte, USA, 2004

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validOptions.insert("prefix", "name");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    interTrackManager manager(args, runTime, mesh);

    manager.read();
    manager.init();

    SM_GLOBALREGIONSCOPE(interTrackManager, DEFAULT);

    volScalarField& p = storage.p();
    volVectorField& U = storage.U();
    surfaceScalarField& phi = storage.phi();
    volScalarField& rho = storage.rho();
    volScalarField& mu = storage.mu();
    trackedSurface& interface = storage.interface();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// TODO: Test F
    storage.item_F().enable();
    volVectorField& F = storage.F();
    F == dimensionedVector
    (
        word(),
        dimForce/dimVolume,
        vector(4, -9.81, 4)
//         vector(0, -9.81, 0)
    );
    Info << F.boundaryField() << endl;

// TODO: Two fluids
// #   include "setRefCell.H"

// TODO: Continuity errors into manager?
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

// TODO: Move to manager function called something like
//       loopSurfaceStatistics
#           include "freeSurfaceContinuityErrs.H"
        }

#       include "volContinuity.H"

// TODO: Move to surface info from writePre to manager
//       function called something like timeSurfaceStatistics
    }

    return(0);
}

// ************************************************************************* //
