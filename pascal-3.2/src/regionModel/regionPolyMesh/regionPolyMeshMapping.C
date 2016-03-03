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

#include "directMappedPatchBase.H"
#include "directMappedPolyPatch.H"
#include "directMappedWallPolyPatch.H"

#include "regionPolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void regionPolyMesh::map
(
    pointField& givenPoints,
    const label& regionI,
    const word& patchName
) const
{
    label regionI0 = regionIndex(polyMesh::defaultRegion);

    const pointField& points = mesh(regionI0).points();

    const labelIOList& map = pointMap(regionI);

    if (patchName == "")
    {
        forAll(givenPoints, pointI)
        {
            givenPoints[pointI] = points[map[pointI]];
        }
    }
    else
    {
        const polyBoundaryMesh& boundaryMesh =
            mesh(regionI).boundaryMesh();

        label patchI = boundaryMesh.findPatchID(patchName);

        if (patchI == -1)
        {
            FatalErrorIn("regionPolyMesh::rmap(...)")
                << "Given patch name " << patchName
                << " for point mapping does not exist in"
                << " region " << regionName(regionI)
                << abort(FatalError);
        }

        const polyPatch& patch = boundaryMesh[patchI];

        labelList mpmap = patch.meshPointMap().toc();

        forAll(mpmap, patchPointI)
        {
            label pointI = mpmap[patchPointI];

            givenPoints[pointI] = points[map[pointI]];
        }
    }
}

void regionPolyMesh::map
(
    pointField& givenPoints,
    const word& regionName,
    const word& patchName
) const
{
    map
    (
        givenPoints,
        regionIndex(regionName),
        patchName
    );
}

tmp<pointField> regionPolyMesh::map
(
    const label& regionI,
    const word& patchName
) const
{
    tmp<pointField> tNewPoints
    (
        new pointField(mesh(regionI).points())
    );

    pointField& newPoints = tNewPoints();

    map
    (
        newPoints,
        regionI,
        patchName
    );

    return tNewPoints;
}

tmp<pointField> regionPolyMesh::map
(
    const word& regionName,
    const word& patchName
) const
{
    return map
    (
        regionIndex(regionName),
        patchName
    );
}

void regionPolyMesh::rmap
(
    pointField& givenPoints,
    const label& regionI,
    const word& patchName
) const
{
    const pointField& points = mesh(regionI).points();

    const labelIOList& map = pointMap(regionI);

    if (patchName == "")
    {
        forAll(points, pointI)
        {
            givenPoints[map[pointI]] = points[pointI];
        }
    }
    else
    {
        const polyBoundaryMesh& boundaryMesh =
            mesh(regionI).boundaryMesh();

        label patchI = boundaryMesh.findPatchID(patchName);

        if (patchI == -1)
        {
            FatalErrorIn("regionPolyMesh::rmap(...)")
                << "Given patch name " << patchName
                << " for point mapping does not exist in"
                << " region " << regionName(regionI)
                << abort(FatalError);
        }

        const polyPatch& patch = boundaryMesh[patchI];

        labelList mpmap = patch.meshPointMap().toc();

        forAll(mpmap, patchPointI)
        {
            label pointI = mpmap[patchPointI];

            givenPoints[map[pointI]] = points[pointI];
        }
    }
}

void regionPolyMesh::rmap
(
    pointField& givenPoints,
    const word& regionName,
    const word& patchName
) const
{
    rmap
    (
        givenPoints,
        regionIndex(regionName),
        patchName
    );
}

tmp<pointField> regionPolyMesh::rmap
(
    const label& regionI,
    const word& patchName
) const
{
    label regionI0 = regionIndex(polyMesh::defaultRegion);

    tmp<pointField> tNewPoints
    (
        new pointField(mesh(regionI0).points())
    );

    pointField& newPoints = tNewPoints();

    rmap
    (
        newPoints,
        regionI,
        patchName
    );

    return tNewPoints;
}

tmp<pointField> regionPolyMesh::rmap
(
    const word& regionName,
    const word& patchName
) const
{
    return rmap
    (
        regionIndex(regionName),
        patchName
    );
}

labelListList regionPolyMesh::patchMapDirectMapped
(
    const label& fromRegionI,
    const label& toRegionI
) const
{
    const polyMesh& fromMesh = mesh(fromRegionI);
    const polyMesh& toMesh = mesh(toRegionI);

    label pmapSize = fromMesh.boundaryMesh().size();

    labelListList pmap
    (
        2,
        labelList
        (
            pmapSize,
            -1
        )
    );

    forAll(fromMesh.boundaryMesh(), fromPatchI)
    {
        if
        (
            (
                fromMesh.boundaryMesh()[fromPatchI].type()
            == directMappedPolyPatch::typeName
            )
            ||
            (
                fromMesh.boundaryMesh()[fromPatchI].type()
            == directMappedWallPolyPatch::typeName
            )
        )
        {
            const directMappedPatchBase& dmpb =
                refCast<const directMappedPatchBase>
                (
                    fromMesh.boundaryMesh()[fromPatchI]
                );

            word toPatchName = dmpb.samplePatch();

            label toPatchI =
                toMesh.boundaryMesh().findPatchID(toPatchName);

            pmap[0][fromPatchI] = toPatchI;
            pmap[1][toPatchI] = fromPatchI;
        }
    }

    return pmap;
}

labelListList regionPolyMesh::patchMapDirectMapped
(
    const word& fromRegionName,
    const word& toRegionName
) const
{
    label fromRegionI = regionIndex(fromRegionName);
    label toRegionI = regionIndex(toRegionName);

    return patchMapDirectMapped
    (
        fromRegionI,
        toRegionI
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

