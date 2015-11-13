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

Class
    regionDynamicFvMesh

Description
    ...

SourceFiles
    regionDynamicFvMesh.C
    regionDynamicFvMeshMapping.C

\*---------------------------------------------------------------------------*/

#ifndef regionDynamicFvMesh_H
#define regionDynamicFvMesh_H

#include "dynamicFvMesh.H"

#include "regionFvMesh.H"

// TODO [High]: Make reference to dynamicFvMesh in access operators constant?

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class regionDynamicFvMesh Declaration
\*---------------------------------------------------------------------------*/

class regionDynamicFvMesh
:
    public regionFvMesh
{
private:

    // Private data

        //- Meshes data
        List<autoPtr<dynamicFvMesh>*> meshesData_;

        //- Meshes access pointers
        List<dynamicFvMesh*> meshes_;

        //- Finite element motion solver type
        boolList isFeMotionSolver_;

        //- Finite volume motion solver type
        boolList isFvMotionSolver_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionDynamicFvMesh(const regionDynamicFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const regionDynamicFvMesh&);

        // Helper functions for constructor

            //- Create new mesh for region with index regionI
            autoPtr<dynamicFvMesh>* newMesh(const label& regionI) const;

public:

    //- Runtime type information
    TypeName("regionDynamicFvMesh");

    // Constructors

        //- Construct from time and list of region names
        regionDynamicFvMesh
        (
            const Time& runTime,
            const wordList& regionNames,
            bool init = true
        );

    // Destructor

        virtual ~regionDynamicFvMesh();

    // Member Functions

        // Access

            //- Return reference to mesh access ist
            const List<dynamicFvMesh*>& meshes() const
            {
                return meshes_;
            };

            //- Return reference to default mesh
            dynamicFvMesh& mesh() const
            {
                return *meshes_[regionIndex(polyMesh::defaultRegion)];
            }

            //- Return reference to mesh
            dynamicFvMesh& mesh(const label& regionI) const
            {
                return *meshes_[regionI];
            }

            //- Return reference to mesh
            dynamicFvMesh& mesh(const word& regionName) const
            {
                return *meshes_[regionIndex(regionName)];
            }

            //- Return refrence to fe motion solver list
            const boolList& isFeMotionSolver() const
            {
                return isFeMotionSolver_;
            }

            //- Return true if motion solver type is fe
            bool isFeMotionSolver(const label& regionI) const
            {
                return isFeMotionSolver_[regionI];
            }

            //- Return true if motion solver type is fe
            bool isFeMotionSolver(const word& regionName) const
            {
                return isFeMotionSolver_[regionIndex(regionName)];
            }

            //- Return refrence to fv motion solver list
            const boolList& isFvMotionSolver() const
            {
                return isFvMotionSolver_;
            }

            //- Return true if motion solver type is fv
            bool isFvMotionSolver(const label& regionI) const
            {
                return isFvMotionSolver_[regionI];
            }

            //- Return true if motion solver type is fv
            bool isFvMotionSolver(const word& regionName) const
            {
                return isFvMotionSolver_[regionIndex(regionName)];
            }

        // Mapping

            //- Map boundary field of mesh velocity between
            //  two neighbour-regions based on directMappedPatch
            void patchMapMeshVelocityDirectMapped
            (
                const label& fromRegionI,
                const label& toRegionI
            ) const;

            //- Map boundary field of mesh velocity between
            //  two neighbour-regions based on directMappedPatch
            void patchMapMeshVelocityDirectMapped
            (
                const word& fromRegionName,
                const word& toRegionName
            ) const;

    // Member operators

        //- Mesh access operator for default region (lookup)
        dynamicFvMesh& operator()() const
        {
            return *meshes_[regionIndex(polyMesh::defaultRegion)];
        }

        //- Mesh access operator by index (direct)
        dynamicFvMesh& operator[](const label& regionI) const
        {
            return *meshes_[regionI];
        }

        //- Mesh access operator by regionName (lookup)
        dynamicFvMesh& operator()(const word& regionName) const
        {
            return *meshes_[regionIndex(regionName)];
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
