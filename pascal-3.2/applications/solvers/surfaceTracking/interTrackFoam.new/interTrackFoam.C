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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "addTrackedSurfacePrefixOption.H"

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    interTrackManager manager(args, runTime, mesh);

    manager.init();

    interTrackControl& control = manager.control();

    interTrackManager::storage& data = manager.data();

    uniformDimensionedVectorField& g = data.g();
    volScalarField& p = data.p();
    volVectorField& U = data.U();
    surfaceScalarField& phi = data.phi();
    volScalarField& rho = data.rho();
    volScalarField& mu = data.mu();
    volVectorField& F = data.F();
    volScalarField& fluidIndicator = data.fluidIndicator();
    twoPhaseMixture& transport = data.transport();
    incompressible::turbulenceModel& turbulence = data.turbulence();
    trackedSurface& interface = data.interface();


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// TODO: Two fluids
// #   include "setRefCell.H"

#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    while (manager.run())
    {
        interface.updateMesh();
        interface.updateDisplacementDirections();

        interface.predictPoints();

        // --- Pressure-velocity PIMPLE corrector loop
        while (control.loop())
        {
            // Update interface bc
            interface.updateBoundaryConditions();

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phi,U);

#           include "interTrackUEqn.H"

             // --- Pressure corrector PISO loop
            while (control.correct())
            {
#               include "interTrackpEqn.H"

                // Make the fluxes relative to the mesh motion
                fvc::makeRelative(phi, U);
            }

            // Make the fluxes absolute
            fvc::makeAbsolute(phi, U);

            interface.correctPoints();

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phi, U);

            if (control.turbCorr())
            {
                transport.correct();
                turbulence.correct();
            }

            // Make the fluxes absolute
            fvc::makeAbsolute(phi, U);

// TODO: Move to manager function called something like
//       loopSurfaceStatistics
#           include "freeSurfaceContinuityErrs.H"
        }

#       include "volContinuity.H"

// TODO: Move to surface info from writePre to manager
//       function called something like timeSurfaceStatistics
        manager.write();
    }

    return(0);
}

// ************************************************************************* //
