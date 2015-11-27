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

void regionPolyMesh::initMeshes(const wordList& regionNames) const
{
    size_ = 1 + regionNames.size();
    regionNames_ = regionNames;

    resizeLists();

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

    setParallelSplitRegions();

    initialized_ = true;
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

wordList regionPolyMesh::readRegionNames() const
{
    wordIOList regionNames
    (
        IOobject
        (
            "regions",
            this->time().constant(),
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    return regionNames;
}

void regionPolyMesh::resizeLists() const
{
    meshesData_.resize(size_, NULL);
    meshes_.resize(size_, NULL);

    cellRegionMap_.resize(size_, NULL);
    pointRegionMap_.resize(size_, NULL);
    faceRegionMap_.resize(size_, NULL);
    cellProcMap_.resize(size_, NULL);
    pointProcMap_.resize(size_, NULL);
    faceProcMap_.resize(size_, NULL);
    cellRegionProcMap_.resize(size_, NULL);
    pointRegionProcMap_.resize(size_, NULL);
    faceRegionProcMap_.resize(size_, NULL);
    cellMap_.resize(size_, NULL);
    pointMap_.resize(size_, NULL);
    faceMap_.resize(size_, NULL);
}

void regionPolyMesh::setParallelSplitRegions() const
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
    bool init
)
:
    objectRegistry(
        IOobject
        (
            "regionPolyMesh",
            runTime.constant(),
            runTime.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    time_(runTime),
    parallel_(runTime.processorCase()),
    size_ (0),
    regionNames_(List<word>(0)),
    parallelSplitRegions_(0),
    meshesData_(List<polyMesh*>(0)),
    meshes_(List<polyMesh*>(0)),
    initialized_(false),
    cellRegionMap_(List<labelIOList*>(0)),
    pointRegionMap_(List<labelIOList*>(0)),
    faceRegionMap_(List<labelIOList*>(0)),
    cellProcMap_(List<labelIOList*>(0)),
    pointProcMap_(List<labelIOList*>(0)),
    faceProcMap_(List<labelIOList*>(0)),
    cellRegionProcMap_(List<labelIOList*>(0)),
    pointRegionProcMap_(List<labelIOList*>(0)),
    faceRegionProcMap_(List<labelIOList*>(0)),
    cellMap_(List<labelIOList*>(0)),
    pointMap_(List<labelIOList*>(0)),
    faceMap_(List<labelIOList*>(0))
{
    if(init) initMeshes(readRegionNames());
}

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
            "regionPolyMesh",
            runTime.constant(),
            runTime.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    time_(runTime),
    parallel_(runTime.processorCase()),
    size_ (1 + regionNames.size()),
    regionNames_(regionNames),
    parallelSplitRegions_(0),
    meshesData_(List<polyMesh*>(0)),
    meshes_(List<polyMesh*>(0)),
    initialized_(false),
    cellRegionMap_(List<labelIOList*>(0)),
    pointRegionMap_(List<labelIOList*>(0)),
    faceRegionMap_(List<labelIOList*>(0)),
    cellProcMap_(List<labelIOList*>(0)),
    pointProcMap_(List<labelIOList*>(0)),
    faceProcMap_(List<labelIOList*>(0)),
    cellRegionProcMap_(List<labelIOList*>(0)),
    pointRegionProcMap_(List<labelIOList*>(0)),
    faceRegionProcMap_(List<labelIOList*>(0)),
    cellMap_(List<labelIOList*>(0)),
    pointMap_(List<labelIOList*>(0)),
    faceMap_(List<labelIOList*>(0))
{
    if(init) initMeshes(regionNames);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label regionPolyMesh::regionIndex() const
{
    return regionIndex(polyMesh::defaultRegion);

    // TODO: If default region not allocated, return -1
}

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

    // TODO: Catch default region here and return -1
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

