/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

#include "regionPolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(regionPolyMesh, 0);

word regionPolyMesh::regionPolyMeshSubDir = "regionPolyMesh";

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

polyMesh* regionPolyMesh::newMesh(const label& regionI) const
{
    return new polyMesh
    (
        IOobject
        (
            regionNames_[regionI],
            time_.timeName(),
            time_,
            IOobject::MUST_READ
        )
    );
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void regionPolyMesh::setParallelSplitRegions()
{
    if (parallel())
    {
        forAll(regionNames(), regionI)
        {
            IOobject faceProcAddObj
            (
                "faceProcAddressing",
                time().constant(),
                polyMesh::meshSubDir,
                mesh(regionI),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            );

            IOobject cellProcAddObj
            (
                "cellProcAddressing",
                time().constant(),
                polyMesh::meshSubDir,
                mesh(regionI),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            );

            IOobject pointProcAddObj
            (
                "pointProcAddressing",
                time().constant(),
                polyMesh::meshSubDir,
                mesh(regionI),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            );

            IOobject boundaryProcAddObj
            (
                "boundaryProcAddressing",
                time().constant(),
                polyMesh::meshSubDir,
                mesh(regionI),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            );

            // If any mesh does not contain processor addressing,
            // assume the mesh was split AFTER decomposition
            if
            (
                !faceProcAddObj.headerOk()
             || !cellProcAddObj.headerOk()
             || !pointProcAddObj.headerOk()
             || !boundaryProcAddObj.headerOk()
            )
            {
                parallelSplitRegions_ = 1;
                break;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

regionPolyMesh::regionPolyMesh
(
    const Time& runTime,
    const wordList& regionNames,
    bool init
)
:
    objectRegistry(
        IOobject
        (
            "multiRegion",
            runTime.constant(),
            runTime.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    time_(runTime),
    size_ (1 + regionNames.size()),
    regionNames_(regionNames),
    parallel_(runTime.processorCase()),
    parallelSplitRegions_(0),
    meshesData_(List<polyMesh*>(size_,NULL)),
    meshes_(List<polyMesh*>(size_,NULL)),
    cellRegionMap_(List<labelIOList*>(size_,NULL)),
    pointRegionMap_(List<labelIOList*>(size_,NULL)),
    faceRegionMap_(List<labelIOList*>(size_,NULL)),
    cellProcMap_(List<labelIOList*>(size_,NULL)),
    pointProcMap_(List<labelIOList*>(size_,NULL)),
    faceProcMap_(List<labelIOList*>(size_,NULL)),
    cellRegionProcMap_(List<labelIOList*>(size_,NULL)),
    pointRegionProcMap_(List<labelIOList*>(size_,NULL)),
    faceRegionProcMap_(List<labelIOList*>(size_,NULL)),
    cellMap_(List<labelIOList*>(size_,NULL)),
    pointMap_(List<labelIOList*>(size_,NULL)),
    faceMap_(List<labelIOList*>(size_,NULL))
{
    if(init)
    {
        forAll(regionNames_, regionI)
        {
            if (debug)
            {
                Info << "regionPolyMesh::regionPolyMesh(...) : "
                    << "Create mesh for region "
                    << regionName(regionI)
                    << endl;
            }

            // Create mesh
            meshesData_[regionI] = newMesh(regionI);

            // Link access pointer
            meshes_[regionI] = meshesData_[regionI];
        }

        // Set parallel split type
        setParallelSplitRegions();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

regionPolyMesh::~regionPolyMesh()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label regionPolyMesh::regionIndex(const word& regionName) const
{
    label regionID = -1;

    forAll(regionNames(), regionI)
    {
        if(regionNames()[regionI] == regionName)
        {
            regionID = regionI;
        }
    }

    if(regionID == -1)
    {
        FatalErrorIn("regionPolyMesh::regionIndex")
            << "Region " << regionName
            << " is not part of this regionPolyMesh!"
            << abort(FatalError);
    }

    return regionID;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

