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

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "trackedSurface.H"
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

    pimpleControl& pimple = manager.pimple();

    interTrackManager::storage& data = manager.data();

    volScalarField& p = data.p();
    volVectorField& U = data.U();
    surfaceScalarField& phi = data.phi();
    volScalarField& rho = data.rho();
    volScalarField& mu = data.mu();
    volScalarField& fluidIndicator = data.fluidIndicator();
    trackedSurface& interface = data.interface();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// TODO

#   include "createSf.H"

// Read interpolators if present
if (interface.twoFluids())
{
    interface.interpolatorAB();
}


#   include "setRefCell.H"

volScalarField AU
(
    IOobject
    (
        "AU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("AU", dimDensity/dimTime, 0)
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    while (manager.run())
    {
        interface.updateMesh();
        interface.updateDisplacementDirections();

        interface.predictPoints();

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            // Update interface bc
            interface.updateBoundaryConditions();

            // Make the fluxes relative
            fvc::makeRelative(phi,U);

#           include "interTrackUEqn.H"

             // --- Pressure corrector PISO loop
            while (pimple.correct())
            {
// TODO
//                 volScalarField rUA = 1.0/UEqn.A();
//                 U = rUA*UEqn.H();
//                 phi = (fvc::interpolate(U) & mesh.Sf());

                volScalarField AU = UEqn.A();
                volVectorField HU = UEqn.H();

                U = HU/AU;

// TODO: Port this from OpenFOAM 3.0.x
                bool ddtPhiCorr = false;
#               include "calcPhi.H"

// #               include "correctPhiAtInterface.H"

#               include "scalePhi.H"

                // Non-orthogonal pressure corrector loop
                while (pimple.correctNonOrthogonal())
                {
                    fvScalarMatrix pEqn
                    (
//                         fvm::laplacian(rUA, p) == fvc::div(phi)
                        fvm::laplacian(1.0/fvc::interpolate(AU), p)
                     == fvc::div(phi)
                    );

#                   include "setReference.H"
//                     pEqn.setReference(pRefCell, pRefValue);

                    pEqn.solve
                    (
                        mesh.solutionDict().solverDict
                        (
                            p.select(pimple.finalInnerIter())
                        )
                    );

                    if (pimple.finalNonOrthogonalIter())
                    {
                        phi -= pEqn.flux();
                    }
                }

#               include "continuityErrs.H"

                // Momentum corrector
                U -= fvc::grad(p)/AU;
//                 U -= rUA*fvc::grad(p);
                U.correctBoundaryConditions();
            }

            interface.correctPoints();

#           include "freeSurfaceContinuityErrs.H"

#           include "updateSf.H"

            AU = UEqn.A();
        }

#       include "volContinuity.H"

        Info << "Total surface tension force: "
            << interface.totalSurfaceTensionForce() << endl;

        vector totalForce =
            interface.totalViscousForce()
          + interface.totalPressureForce();

        Info << "Total force: " << totalForce << endl;

// TEST: Sub-mesh
        const areaScalarField& K = interface.curvature();
//        const areaScalarField& K = interface.aMesh().faceCurvatures();

        Info << "Free surface curvature: min = " << gMin(K)
            << ", max = " << gMax(K)
            << ", average = " << gAverage(K) << endl << flush;

        manager.write();
    }

    return(0);
}

// ************************************************************************* //
