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

#include "tetPointFields.H"
#include "pointFields.H"

#include "regionDynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(regionDynamicFvMesh, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

autoPtr<dynamicFvMesh>*
regionDynamicFvMesh::newMesh(const label& regionI) const
{
    return new autoPtr<dynamicFvMesh>
    (
        dynamicFvMesh::New
        (
            IOobject
            (
                regionNames_[regionI],
                time_.timeName(),
                time_,
                IOobject::MUST_READ
            )
        )
    );
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void regionDynamicFvMesh::resizeLists() const
{
    regionFvMesh::resizeLists();

    isFeMotionSolver_.resize(size_, NULL);
    isFvMotionSolver_.resize(size_, NULL);

    dynamicFvMeshesData_.resize(size_, NULL);
    dynamicFvMeshes_.resize(size_, NULL);
}

void regionDynamicFvMesh::initMeshes(const wordList& regionNames) const
{
    size_ = 1 + regionNames.size();
    regionNames_ = regionNames;

    resizeLists();

    forAll(regionNames_, regionI)
    {
        if (debug)
        {
            Info << "regionDynamicFvMesh::regionDynamicFvMesh(...) : "
                << "Create mesh for region "
                << regionName(regionI)
                << endl;
        }

        // Create mesh
        dynamicFvMeshesData_[regionI] = newMesh(regionI);

        // Link access pointer
        dynamicFvMeshes_[regionI] = dynamicFvMeshesData_[regionI]->operator->();

        // Remember if type of motion solver is fe
        isFeMotionSolver_[regionI] = dynamicFvMeshes_[regionI]
            ->objectRegistry::foundObject<tetPointVectorField>
            (
                "motionU"
            );

        // Remember if type of motion solver is fv
        isFvMotionSolver_[regionI] = dynamicFvMeshes_[regionI]
            ->objectRegistry::foundObject<pointVectorField>
            (
                "pointMotionU"
            );
    }

    setParallelSplitRegions();

    initialized_ = true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

regionDynamicFvMesh::regionDynamicFvMesh
(
    const Time& runTime,
    bool init
)
:
    regionFvMesh::regionFvMesh
    (
        runTime,
        false
    ),
    dynamicFvMeshesData_(List<autoPtr<dynamicFvMesh>*>(0)),
    dynamicFvMeshes_(List<dynamicFvMesh*>(0)),
    isFeMotionSolver_(boolList(0)),
    isFvMotionSolver_(boolList(0))
{
    if(init) initMeshes(readRegionNames());
}

regionDynamicFvMesh::regionDynamicFvMesh
(
    const Time& runTime,
    const wordList& regionNames,
    bool init
)
:
    regionFvMesh::regionFvMesh
    (
        runTime,
        false
    ),
    dynamicFvMeshesData_(List<autoPtr<dynamicFvMesh>*>(0)),
    dynamicFvMeshes_(List<dynamicFvMesh*>(0)),
    isFeMotionSolver_(boolList(0)),
    isFvMotionSolver_(boolList(0))
{
    if(init) initMeshes(regionNames);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

