/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Application
    interfaceTrackinFoam

Description
    Incompressible laminar CFD code for simulation of a single bubble rising 
    in a stil liquid. Interface between fluid phases is tracked using moving
    mesh.

\*---------------------------------------------------------------------------*/

#include "Ostream.H"
#include "OFstream.H"

#include "fvCFD.H"

#include "zeroGradientFvPatchFields.H"

#include "freeSurface.H"
#include "dynamicFvMesh.H"

#include "tetFemMatrices.H"
#include "tetPointFields.H"
#include "faceTetPolyPatch.H"
#include "tetPolyPatchInterpolation.H"
#include "fixedValueTetPolyPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"

#   include "createDynamicFvMesh.H"

#   include "createFields.H"

    // Spherical bubble radius
//     scalar r = 0.001;

    // Oblate spheroid dimensions
    scalar a = 0.00105;
    scalar b = 0.0009759;
    scalar c = b;

    freeSurface interface(mesh, rho, U, p, phi);

    runTime++;

    Info << "Adjust interface points ...";

    // Check mesh motion solver type 
    bool feMotionSolver = 
        mesh.objectRegistry::foundObject<tetPointVectorField>
        (
            "motionU"
        );
//     bool fvMotionSolver =
//         mesh.objectRegistry::foundObject<pointVectorField>
//         (
//             "pointMotionU"
//         );

    pointField newPoints = mesh.points();

    if (feMotionSolver)
    {
        tetPointVectorField& motionU =
            const_cast<tetPointVectorField&>
            (
                mesh.objectRegistry::
                lookupObject<tetPointVectorField>
                (
                    "motionU"
                )
            );

        fixedValueTetPolyPatchVectorField& motionUaPatch =
            refCast<fixedValueTetPolyPatchVectorField>
            (
                motionU.boundaryField()[interface.aPatchID()]
            );

        tetPolyPatchInterpolation tppiAPatch
        (
            refCast<const faceTetPolyPatch>
            (
                motionUaPatch.patch()
            )
        );

        const labelList& meshPointsA = 
            mesh.boundaryMesh()[interface.aPatchID()].meshPoints();

        vectorField totalDisplacement(meshPointsA.size(), vector::zero);

        forAll(meshPointsA, pointI)
        {
            label curPoint = meshPointsA[pointI];

            scalar sqrSinV = 
                sqr(newPoints[curPoint].z()/mag(newPoints[curPoint]));

            scalar sqrCosUCosV =
                sqr(newPoints[curPoint].x()/mag(newPoints[curPoint]));

            scalar sqrCosVSinU =
                sqr(newPoints[curPoint].y()/mag(newPoints[curPoint]));

            scalar newR = 
                1.0
               /::sqrt
                (
                    sqrCosUCosV/sqr(a)
                  + sqrCosVSinU/sqr(b)
                  + sqrSinV/sqr(c)
                );

            totalDisplacement[pointI] =
                (newR/mag(newPoints[curPoint]) - 1.0)
               *newPoints[curPoint];
        }    

        motionUaPatch ==
            tppiAPatch.pointToPointInterpolate
            (
                totalDisplacement/runTime.deltaT().value()
            );

        if(interface.twoFluids())
        {
            pointField totDisplacementB =
                interface.interpolatorAB().pointInterpolate
                (
                    totalDisplacement
                );
        
            fixedValueTetPolyPatchVectorField& motionUbPatch =
                refCast<fixedValueTetPolyPatchVectorField>
                (
                    motionU.boundaryField()[interface.bPatchID()]
                );

            tetPolyPatchInterpolation tppiBPatch
                (
                    refCast<const faceTetPolyPatch>(motionUbPatch.patch())
                );

            motionUbPatch == 
                tppiBPatch.pointToPointInterpolate
                (
                    totDisplacementB/runTime.deltaT().value()
                );
        }
    }

    mesh.update();

    interface.aMesh().movePoints(mesh.points());

    Info << "done" << endl;

    if (interface.twoFluids())
    {
        Info << "Create interface interpolators ..." << flush;

        interface.interpolatorAB();

        Info << "done" << endl;
    }

    runTime.writeNow();

    Info << "ExecutionTime = "
        << scalar(runTime.elapsedCpuTime())
        << " s\n" << endl << endl;
    
    Info << "End\n" << endl;

    return(0);
}

// ************************************************************************* //
