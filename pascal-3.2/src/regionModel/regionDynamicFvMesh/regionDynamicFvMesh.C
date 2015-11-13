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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

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
        regionNames,
        false
    ),
    meshesData_(List<autoPtr<dynamicFvMesh>*>(size_,NULL)),
    meshes_(List<dynamicFvMesh*>(size_,NULL)),
    isFeMotionSolver_(boolList(size_,false)),
    isFvMotionSolver_(boolList(size_,false))
{
    if(init)
    {
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
            meshesData_[regionI] = newMesh(regionI);

            // Link access pointer
            meshes_[regionI] = meshesData_[regionI]->operator->();

            // Link access pointer of regionPolyMesh class
            regionPolyMesh::meshes_[regionI] =
                meshes_[regionI];

            // Link access pointer of regionFvMesh class
            regionFvMesh::meshes_[regionI] =
                meshes_[regionI];

            // Remember if type of motion solver is fe
            isFeMotionSolver_[regionI] = meshes_[regionI]
                ->objectRegistry::foundObject<tetPointVectorField>
                (
                    "motionU"
                );

            // Remember if type of motion solver is fv
            isFvMotionSolver_[regionI] = meshes_[regionI]
                ->objectRegistry::foundObject<pointVectorField>
                (
                    "pointMotionU"
                );
        }

        // Set parallel split type
        setParallelSplitRegions();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

regionDynamicFvMesh::~regionDynamicFvMesh()
{
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

