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
    solverManagerRegions

Description
    ...

\*---------------------------------------------------------------------------*/

#ifndef solverManagerRegions_H
#define solverManagerRegions_H

#include "autoPtr.H"
#include "argList.H"
#include "foamTime.H"
#include "IOdictionary.H"
#include "solutionControl.H"
#include "dynamicFvMesh.H"

#include "solverManager.H"
#include "solverManagerRegion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define solverManagerRegions_new(MANAGER, NAME)                              \
                                                                              \
    solverManagerNestedClass_new(MANAGER, Regions, NAME)                      \
                                                                              \
    protected: virtual void init(const word init = "default") const {};


# define solverManagerRegions_addDefaultRegion(TYPE)                          \
                                                                              \
    solverManagerNestedClass_addVirtual(TYPE, defaultRegion)


# define solverManagerRegions_addRegion(TYPE, REGION)                         \
                                                                              \
    solverManagerNestedClass_addVirtual(TYPE, REGION)


# define solverManagerRegions_end()                                           \
                                                                              \
    solverManagerNestedClass_end()


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class solverManagerRegions Declaration
\*---------------------------------------------------------------------------*/

template <class MESH>
class solverManagerRegions
{
public:

    // Public typedefs

        typedef MESH Mesh;

        typedef solverManagerRegion<MESH> Region;
        typedef solutionControl Control;
        typedef solverManagerStorage<MESH> Storage;


private:

    // Private Data

        //- Reference to args
        const argList& args_;

        //- Reference to time
        const Time& time_;

        //- Reference to mesh
        const MESH& mesh_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        solverManagerRegions(const solverManagerRegions&);

        //- Disallow default bitwise assignment
        void operator=(const solverManagerRegions&);


protected:

    // Protected Member Functions

        //- Reference to args
        inline const argList& args() const
        {
            return args_;
        };

        //- Reference to time
        inline Time& time() const
        {
            return time_;
        };

        //- Return reference to mesh
        inline MESH& mesh() const
        {
            return mesh_;
        }

        //- Initialize data
        virtual void init(const word init) const = 0;


public:

    // Constructor

        solverManagerRegions
        (
            const argList& args,
            const Time& time,
            const MESH& mesh
        )
        :
            args_(args), time_(time), mesh_(mesh)
        {};

    // Destructor

        virtual ~solverManagerRegions() {};

    // Member Functions

        //- Default region
        virtual Region& defaultRegion() const = 0;
};


template <>
class solverManagerRegions<dynamicFvMesh>
{
public:

    // Public typedefs

        typedef dynamicFvMesh Mesh;

        typedef solverManagerRegion<dynamicFvMesh> Region;
        typedef solutionControl Control;
        typedef solverManagerStorage<dynamicFvMesh> Storage;


private:

    // Private Data

        //- Reference to args
        const argList& args_;

        //- Reference to time
        Time& time_;

        //- Reference to mesh
        dynamicFvMesh& mesh_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        solverManagerRegions(const solverManagerRegions&);

        //- Disallow default bitwise assignment
        void operator=(const solverManagerRegions&);


protected:

    // Protected Member Functions

        //- Reference to args
        inline const argList& args() const
        {
            return args_;
        };

        //- Reference to time
        inline Time& time() const
        {
            return time_;
        };

        //- Return reference to mesh
        inline dynamicFvMesh& mesh() const
        {
            return mesh_;
        }

        //- Initialize data
        virtual void init(const word init) const = 0;


public:

    // Constructor

        solverManagerRegions
        (
            const argList& args,
            Time& time,
            dynamicFvMesh& mesh
        )
        :
            args_(args), time_(time), mesh_(mesh)
        {};

    // Destructor

        virtual ~solverManagerRegions() {};

    // Member Functions

        //- Default region
        virtual Region& defaultRegion() const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
