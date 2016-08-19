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
    regionPointMesh

Description
    ...

SourceFiles
    regionPointMesh.C

\*---------------------------------------------------------------------------*/

#ifndef regionPointMesh_H
#define regionPointMesh_H

#include "pointMesh.H"

#include "regionPolyMesh.H"

// TODO [Low]: Implement mapping

// TODO [Low]: Derive regionPointMesh from regionGeoMesh<regionPolyMesh>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class regionPointMesh Declaration
\*---------------------------------------------------------------------------*/

class regionPointMesh
{
private:

    // Private data

        //- Reference to regionPolyMesh
        const regionPolyMesh& rpMesh_;

        //- Mesh pointers
        mutable PtrList<pointMesh> meshPtrs_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionPointMesh(const regionPointMesh&);

        //- Disallow default bitwise assignment
        void operator=(const regionPointMesh&);

        // Helper functions for constructor

            //- Create new mesh for region with index regionI
            pointMesh* newMesh(label regionI) const;


public:

    // Public typedefs

        typedef pointMesh Mesh;
        typedef pointBoundaryMesh BoundaryMesh;

    //- Runtime type information
    TypeName("regionPointMesh");

    // Constructors

        //- Construct from regionPolyMesh
        explicit regionPointMesh
        (
            const regionPolyMesh& rpMesh
        );

    // Destructor

        virtual ~regionPointMesh();

    // Member Functions

        // Access

            //- Return reference to regionPolyMesh
            inline const regionPolyMesh& rpmesh() const
            {
                return rpMesh_;
            };

            //- Return reference to this mesh
            inline const regionPointMesh& rmesh() const
            {
                return *this;
            }

            //- Return the object registry
            inline const objectRegistry& thisDb() const
            {
                return rpMesh_.thisDb();
            }

            //- Return const reference to time
            inline const Time& time() const
            {
                return rpMesh_.time();
            };

            //- Return const reference to global region count
            inline label size() const
            {
                return rpMesh_.size();
            };

            //- Return const reference to regions
            inline const regionList& regions() const
            {
                return rpMesh_.regions();
            };

            //- Return reference to mesh
            inline const pointMesh& mesh
            (
                label regionI
            ) const
            {
                return operator[](regionI);
            }

            //- Return reference to mesh
            inline const pointMesh& mesh
            (
                const word& regionName = polyMesh::defaultRegion
            ) const
            {
                return operator[](regions()[regionName]);
            }

    // Member Operators

        //- Mesh access operator by index (direct)
        inline pointMesh& operator[]
        (
            label regionI
        ) const
        {
            return meshPtrs_[regionI];
        }

        //- Mesh access operator by regionName
        inline pointMesh& operator[]
        (
            const word& regionName
        ) const
        {
            return operator[](regions()[regionName]);
        }

        bool operator!=(const regionPointMesh&) const;
        bool operator==(const regionPointMesh&) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
