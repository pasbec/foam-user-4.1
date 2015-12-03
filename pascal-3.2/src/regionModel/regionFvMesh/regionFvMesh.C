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

#include "regionFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(regionFvMesh, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

fvMesh* regionFvMesh::newMesh(const label& regionI) const
{
    return new fvMesh
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


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void regionFvMesh::initMeshes(const wordList& regionNames) const
{
    size_ = 1 + regionNames.size();
    regionNames_ = regionNames;

    resizeLists();

    forAll(regionNames_, regionI)
    {
        if (debug)
        {
            Info << "regionFvMesh::regionFvMesh(...) : "
                << "Create mesh for region "
                << regionName(regionI)
                << endl;
        }

        // Create mesh
        fvMeshes_[regionI] = newMesh(regionI);
    }

    setParallelSplitRegions();

    initialized_ = true;
}

void regionFvMesh::resizeLists() const
{
    regionPolyMesh::resizeLists();

    fvMeshes_.resize(size_, NULL);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

regionFvMesh::regionFvMesh
(
    const Time& runTime,
    bool init
)
:
    regionPolyMesh::regionPolyMesh
    (
        runTime,
        false
    ),
    fvMeshes_(List<fvMesh*>(0))
{
    if(init) initMeshes(readRegionNames());
}

regionFvMesh::regionFvMesh
(
    const Time& runTime,
    const wordList& regionNames,
    bool init
)
:
    regionPolyMesh::regionPolyMesh
    (
        runTime,
        false
    ),
    fvMeshes_(List<fvMesh*>(0))
{
    if(init) initMeshes(regionNames);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool regionFvMesh::operator!=(const regionFvMesh& brm) const
{
    return &brm != this;
}

bool regionFvMesh::operator==(const regionFvMesh& brm) const
{
    return &brm == this;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

