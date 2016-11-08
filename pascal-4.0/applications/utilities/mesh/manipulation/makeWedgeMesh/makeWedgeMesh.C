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
    makeWedgeMesh

Description
    Make wedge mesh from a 2D mesh

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::validOptions.insert("overwrite", "");
    argList::validOptions.insert("angle", "scalar [deg]");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    bool overwrite = args.optionFound("overwrite");
    word pointsInstance = mesh.pointsInstance();

    scalar angle = 5.0;
    args.optionReadIfPresent("angle", angle);

    scalar phi = mathematicalConstant::pi/180.0 * angle/2.0;

    label nGeometricD = mesh.nGeometricD();
    Vector<label> geometricD = mesh.geometricD();

    labelHashSet nSet;
    labelHashSet tSet;

    if (nGeometricD == 2)
    {
        for (direction cmpt = 0; cmpt < Vector<label>::nComponents; cmpt++)
        {
            if (geometricD[cmpt] < 0)
            {
                nSet.insert(cmpt);
            }
            else
            {
                tSet.insert(cmpt);
            }
        }

        labelList n = nSet.sortedToc();
        labelList t = tSet.sortedToc();

        // Create the new points
        pointField newPoints = mesh.allPoints();

        Info << "Transform points" << endl;

        const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();

        bool mesh2D = false;

        forAll (boundaryMesh, patchI)
        {
            const polyPatch& patch = boundaryMesh[patchI];

            if (isA<emptyPolyPatch>(patch))
            {
                const labelList& meshPoints = patch.meshPoints();

                forAll (meshPoints, patchPointI)
                {
                    label pointI = meshPoints[patchPointI];

                    point p = newPoints[pointI];

                    scalar s = Foam::sign(p[n[0]]);
                    scalar r = p[t[0]];

                    p[n[0]] = s*r*Foam::sin(phi);
                    p[t[0]] =   r*Foam::cos(phi);

                    newPoints[pointI] = p;
                }

                mesh2D = true;

                // This is a cruel hack in order to avoid juggling with
                // the boundaryMesh or the boundary dictionary itself
                const_cast<word&>(patch.type()) = "wedge";
            }
        }

        if (mesh2D)
        {
            if (!overwrite)
            {
                runTime++;

                mesh.setInstance(mesh.time().timeName());
            }

            // This is a cruel hack in order to avoid juggling with
            // the boundaryMesh or the boundary dictionary itself
            mesh.boundaryMesh().write();


            Info << "Write new mesh" << endl;

            // Write mesh (points)
            mesh.write();

//             polyMesh wedgeMesh
//             (
//                 IOobject
//                 (
//                     mesh.name(),
//                     mesh.pointsInstance(),
//                     runTime,
//                     IOobject::NO_READ,
//                     IOobject::NO_WRITE
//                 ),
//                 xferCopy(mesh.allPoints()),
//                 xferCopy(mesh.allFaces()),
//                 xferCopy(mesh.faceOwner()),
//                 xferCopy(mesh.faceNeighbour()),
//                 false
//             );
//
//             Info << "Create new patches" << endl;
//
//             // Create new patches
//             List<polyPatch*> newPatches(boundaryMesh.size(), NULL);
//
//             forAll (boundaryMesh, patchI)
//             {
//                 const polyPatch& patch = boundaryMesh[patchI];
//
//                 word curPatchType = patch.type();
//
//                 // Replace empty patches with wedge
//                 if (curPatchType == emptyPolyPatch::typeName)
//                 {
//                     curPatchType = wedgePolyPatch::typeName;
//                 }
//
//     Info << "patch.name() = " << patch.name() << endl;
//     Info << "patch.type() = " << patch.type() << endl;
//     Info << "patch.start() = " << patch.start() << endl;
//     Info << "patch.size() = " << patch.size() << endl;
//     Info << "curPatchType = " << curPatchType << endl;
//
//                 newPatches[patchI] = polyPatch::New
//                 (
//                     curPatchType,
//                     patch.name(),
//                     patch.size(),
//                     patch.start(),
//                     patchI,
//                     wedgeMesh.boundaryMesh()
//                 ).ptr();
//             }
//
//             // Add new patches
//             wedgeMesh.addPatches(newPatches);
//
//             if (!overwrite)
//             {
//                 runTime++;
//
//                 wedgeMesh.setInstance(mesh.time().timeName());
//             }
//
//             Info << "Write new mesh" << endl;
//
//             wedgeMesh.write();
        }
        else
        {
            Info << "Mesh is 2D, but does not have empty patches!" << endl;
        }
    }
    else
    {
        Info << "Mesh is not 2D!" << endl;
    }

    Info << "\nEnd\n" << endl;

    return(0);
}

// ************************************************************************* //
