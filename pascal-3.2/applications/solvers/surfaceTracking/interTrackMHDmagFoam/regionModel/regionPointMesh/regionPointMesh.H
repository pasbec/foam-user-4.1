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

// TODO [High]: Derive regionPointMesh from regionGeoMesh<regionPolyMesh>

// TODO [High]: Implement access functions for point mapping

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

        //- Meshes
        List<pointMesh*> meshes_;

    // Private demand-driven data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionPointMesh(const regionPointMesh&);

        //- Disallow default bitwise assignment
        void operator=(const regionPointMesh&);

        // Helper functions for constructor

            //- Create new mesh for region with index regionI
            pointMesh* newMesh(const label& regionI) const;

public:

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
            const regionPolyMesh& rpmesh() const
            {
                return rpMesh_;
            };

            //- Return const reference to time
            const Time& time() const
            {
                return rpMesh_.time();
            };

            //- Return const reference to global region count
            const label& size() const
            {
                return rpMesh_.size();
            };

            //- Return const reference to the list of all region names
            const List<word>& regionNames() const
            {
                return rpMesh_.regionNames();
            };

            //- Return name of region with index regionI
            word regionName(const label& regionI) const
            {
                return rpMesh_.regionName(regionI);
            };

            //- Return region index for region with regionName
            label regionIndex(const word& regionName) const
            {
                return rpMesh_.regionIndex(regionName);
            };

            //- Return reference to mesh access list
            const List<pointMesh*>& meshes() const
            {
                return meshes_;
            };

            //- Return reference to default mesh
            const pointMesh& mesh() const
            {
                return *meshes_[regionIndex(polyMesh::defaultRegion)];
            }

            //- Return reference to mesh
            const pointMesh& mesh(const label& regionI) const
            {
                return *meshes_[regionI];
            }

            //- Return reference to mesh
            const pointMesh& mesh(const word& regionName) const
            {
                return *meshes_[regionIndex(regionName)];
            }

            //- Return the object registry
            const objectRegistry& thisDb() const
            {
                return rpMesh_.thisDb();
            }

        // Mapping

        // Demand-driven access

            //- Return const reference to mesh cell mapping
            const labelIOList& cellMap(const label& regionI) const
            {
                return rpMesh_.cellMap(regionI);
            };

            //- Return const reference to mesh cell mapping
            const labelIOList& cellMap(const word& regionName) const
            {
                return rpMesh_.cellMap(regionName);
            };

            //- Return const reference to mesh point mapping
            const labelIOList& pointMap(const label& regionI) const
            {
                return rpMesh_.pointMap(regionI);
            };

            //- Return const reference to mesh point mapping
            const labelIOList& pointMap(const word& regionName) const
            {
                return rpMesh_.pointMap(regionName);
            };

        // Operators

            //- Mesh access operator for default region (lookup)
            const pointMesh& operator()() const
            {
                return *meshes_[regionIndex(polyMesh::defaultRegion)];
            }

            //- Mesh access operator by index (direct)
            const pointMesh& operator[](const label& regionI) const
            {
                return *meshes_[regionI];
            }

            //- Mesh access operator by regionName (lookup)
            const pointMesh& operator()(const word& regionName) const
            {
                return *meshes_[regionIndex(regionName)];
            }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
