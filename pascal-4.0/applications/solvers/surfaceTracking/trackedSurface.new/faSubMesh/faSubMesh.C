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

\*---------------------------------------------------------------------------*/

#include "faSubMesh.H"

#include "triSurface.H"
#include "Xfer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(faSubMesh, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

tmp<pointField> faSubMesh::calcNewPoints() const
{
    const pointField& basePolyMeshPoints = basePolyMesh().points();
    const pointField& baseAreaMeshPoints = baseAreaMesh().points();
    const pointField& splitPoints = subFacePoints();

    tmp<pointField> tpoints
    (
        new pointField(subPolyMesh().allPoints())
    );
    pointField& points = tpoints();

// TODO: Optimize without if statements by knowing the sizes of each sublist
    forAll(subPolyMesh().points(), pointI)
    {
        label basePolyMeshPointI = subPolyMeshPointMap_[pointI];
        label baseAreaMeshPointI = subPolyMeshPointAreaMap_[pointI];
        label splitPointI = subPolyMeshPointTriMap_[pointI];

        if (basePolyMeshPointI > -1)
        {
            points[pointI] = basePolyMeshPoints[basePolyMeshPointI];
        }
        else if (baseAreaMeshPointI > -1)
        {
            points[pointI] = baseAreaMeshPoints[baseAreaMeshPointI];
        }
        else if (splitPointI > -1)
        {
            points[pointI] = splitPoints[splitPointI];
        }
        else
        {
            // FatalError
        }
    }

    return tpoints;
}


void faSubMesh::clearGeom() const
{
    deleteDemandDrivenData(faceCurvaturesPtr_);
}


void faSubMesh::clearOut() const
{
    clearGeom();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

faSubMesh::faSubMesh
(
    const faMesh& baseAreaMesh,
    const pointField& subFacePoints
)
:
basePolyMesh_(baseAreaMesh.mesh()),
baseAreaMesh_(baseAreaMesh),
subPolyMesh_
(
    IOobject
    (
        basePolyMesh_.name() + "_faSubMesh",
        basePolyMesh_.time().timeName(),
        basePolyMesh_.time(),
        IOobject::MUST_READ
    )
),
subPolyMeshPointMap_
(
    IOobject
    (
        "pointMap",
        subPolyMesh_.facesInstance(),
        polyMesh::meshSubDir,
        subPolyMesh_,
        IOobject::MUST_READ
    )
),
subPolyMeshPointAreaMap_
(
    IOobject
    (
        "pointAreaMap",
        subPolyMesh_.facesInstance(),
        polyMesh::meshSubDir,
        subPolyMesh_,
        IOobject::MUST_READ
    )
),
subPolyMeshPointTriMap_
(
    IOobject
    (
        "pointTriMap",
        subPolyMesh_.facesInstance(),
        polyMesh::meshSubDir,
        subPolyMesh_,
        IOobject::MUST_READ
    )
),
subPolyMeshFaceMap_
(
    IOobject
    (
        "faceMap",
        subPolyMesh_.facesInstance(),
        polyMesh::meshSubDir,
        subPolyMesh_,
        IOobject::MUST_READ
    )
),
subPolyMeshCellMap_
(
    IOobject
    (
        "cellMap",
        subPolyMesh_.facesInstance(),
        polyMesh::meshSubDir,
        subPolyMesh_,
        IOobject::MUST_READ
    )
),
subAreaMesh_(subPolyMesh_),
subFacePoints_(subFacePoints),
faceCurvaturesPtr_(NULL)
{}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

faSubMesh::~faSubMesh()
{
    clearOut();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
