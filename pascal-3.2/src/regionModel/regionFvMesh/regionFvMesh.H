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
    regionFvMesh

Description
    ...

SourceFiles
    regionFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef regionFvMesh_H
#define regionFvMesh_H

#include "fvMesh.H"

#include "regionPolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class regionFvMesh Declaration
\*---------------------------------------------------------------------------*/

class regionFvMesh
:
    public regionPolyMesh
{
private:

    // Private data

        //- Meshes data
        mutable List<fvMesh*> fvMeshes_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionFvMesh(const regionFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const regionFvMesh&);

        // Helper functions for constructor

            //- Create new mesh for region with index regionI
            fvMesh* newMesh(const label& regionI) const;

protected:

    // Protected data

    // Private Member Functions

        // Helper functions for constructor

            //- Resize lists storing data
            virtual void resizeLists() const;

            //- Init regions
            virtual void initMeshes(const wordList& regionNames) const;

public:

    //- Runtime type information
    TypeName("regionFvMesh");

    // Constructors

        //- Construct from time
        regionFvMesh
        (
            const Time& runTime,
            bool init = true
        );

        //- Construct from time and list of region names
        regionFvMesh
        (
            const Time& runTime,
            const wordList& regionNames,
            bool init = true
        );

    // Destructor

        virtual ~regionFvMesh() {};

    // Member Functions

        // Access

            //- Return reference to mesh
            fvMesh& mesh(const label& regionI) const
            {
                return operator[](regionI);
            }

            //- Return reference to mesh
            fvMesh& mesh(const word& regionName) const
            {
                return operator[](regionIndex(regionName));
            }

            //- Return reference to default solution dict
            fvSolution& solutionDict() const
            {
                return mesh(regionIndex(polyMesh::defaultRegion)).solutionDict();
            };

            //- Return reference to solution dict
            fvSolution& solutionDict(const label& regionI) const
            {
                return mesh(regionI).solutionDict();
            };

            //- Return reference to solution dict
            fvSolution& solutionDict(const word& regionName) const
            {
                return mesh(regionIndex(regionName)).solutionDict();
            };

    // Member operators

        //- Mesh access operator by index
        virtual fvMesh& operator[](const label& regionI) const
        {
            return *fvMeshes_[regionI];
        }

        //- Mesh access operator by regionName
        polyMesh& operator[](const word& regionName) const
        {
            return operator[](regionIndex(regionName));
        }

        bool operator!=(const regionFvMesh&) const;
        bool operator==(const regionFvMesh&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
