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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"

#   include "createDynamicFvMesh.H"

#   include "createFields.H"

    scalar r = 0.001;

    freeSurface interface(mesh, rho, U, p, phi);

    volScalarField fluidIndicator
     (
         IOobject
         (
             "fluidIndicator",
             runTime.timeName(),
             mesh,
             IOobject::NO_READ,
             IOobject::NO_WRITE
         ),
         interface.fluidIndicator()
     );

    runTime++;

    Info << "Adjust interface points ...";

    vectorField newPoints = mesh.points();

    const labelList& meshPointsA = 
        mesh.boundaryMesh()[interface.aPatchID()].meshPoints();
    
    forAll(meshPointsA, pointI)
    {
        label curPoint = meshPointsA[pointI];

        newPoints[curPoint] *= r/mag(newPoints[curPoint]);
    }    


    if (interface.twoFluids())
    {
        const labelList& meshPointsB = 
            mesh.boundaryMesh()[interface.bPatchID()].meshPoints();

        forAll(meshPointsB, pointI)
        {
            label curPoint = meshPointsB[pointI];
            
            newPoints[curPoint] *= r/mag(newPoints[curPoint]);
        }
    }

    mesh.movePoints(newPoints);

    interface.aMesh().movePoints();

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
