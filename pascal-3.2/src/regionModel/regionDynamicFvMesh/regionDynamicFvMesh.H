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

// TODO [Low]: What if neighbour regions do se different motion solver types?

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
        mutable List<autoPtr<dynamicFvMesh>*> dynamicFvMeshesData_;

        //- Meshes access pointers
        mutable List<dynamicFvMesh*> dynamicFvMeshes_;

        //- Finite element motion solver type
        mutable boolList isFeMotionSolver_;

        //- Finite volume motion solver type
        mutable boolList isFvMotionSolver_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionDynamicFvMesh(const regionDynamicFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const regionDynamicFvMesh&);

        // Helper functions for constructor

            //- Create new mesh for region with index regionI
            autoPtr<dynamicFvMesh>* newMesh(const label& regionI) const;


protected:

    // Protected Member Functions

        // Helper functions for constructor

            //- Resize lists storing data
            virtual void resizeLists() const;

            //- Init regions
            virtual void initMeshes(const wordList& regionNames) const;


public:

    //- Runtime type information
    TypeName("regionDynamicFvMesh");

    // Constructors

        //- Construct from time (without init)
        regionDynamicFvMesh
        (
            const Time& runTime,
            bool init = true
        );

        //- Construct from time and list of region names
        regionDynamicFvMesh
        (
            const Time& runTime,
            const wordList& regionNames,
            bool init = true
        );

    // Destructor

        virtual ~regionDynamicFvMesh() {};

    // Member Functions

        // Access

            //- Return reference to mesh
            inline dynamicFvMesh& mesh(const label& regionI) const
            {
                return operator[](regionI);
            }

            //- Return reference to mesh
            inline dynamicFvMesh& mesh(const word& regionName) const
            {
                return operator[](regionIndex(regionName));
            }

            //- Return refrence to fe motion solver list
            inline const boolList& isFeMotionSolver() const
            {
                return isFeMotionSolver_;
            }

            //- Return true if motion solver type is fe
            inline bool isFeMotionSolver(const label& regionI) const
            {
                return isFeMotionSolver_[regionI];
            }

            //- Return true if motion solver type is fe
            inline bool isFeMotionSolver(const word& regionName) const
            {
                return isFeMotionSolver_[regionIndex(regionName)];
            }

            //- Return refrence to fv motion solver list
            inline const boolList& isFvMotionSolver() const
            {
                return isFvMotionSolver_;
            }

            //- Return true if motion solver type is fv
            inline bool isFvMotionSolver(const label& regionI) const
            {
                return isFvMotionSolver_[regionI];
            }

            //- Return true if motion solver type is fv
            inline bool isFvMotionSolver(const word& regionName) const
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

        //- Mesh access operator by index
        virtual dynamicFvMesh& operator[](const label& regionI) const
        {
            return *dynamicFvMeshes_[regionI];
        }

        //- Mesh access operator by regionName
        inline dynamicFvMesh& operator[](const word& regionName) const
        {
            return operator[](regionIndex(regionName));
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
