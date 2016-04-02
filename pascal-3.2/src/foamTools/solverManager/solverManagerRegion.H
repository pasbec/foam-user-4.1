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
    solverManagerRegion

Description
    ...

\*---------------------------------------------------------------------------*/

#ifndef solverManagerRegion_H
#define solverManagerRegion_H

#include "autoPtr.H"
#include "argList.H"
#include "foamTime.H"
#include "IOdictionary.H"
#include "solutionControl.H"

#include "solverManager.H"
#include "solverManagerStorage.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// TODO: Macro new for regionMeshes where mesh is the mesh of one region

# define solverManagerRegion_new(MANAGER, NAME)                               \
                                                                              \
    solverManagerNestedClass_new(MANAGER, Region, NAME)                       \
                                                                              \
    protected: virtual void init(const word& init = "default") const {};


// TODO: Constructor for control only takes one parameter: mesh
# define solverManagerRegion_addControl(TYPE)                                 \
                                                                              \
    public: typedef TYPE Control;                                             \
                                                                              \
    private: mutable autoPtr<Control> controlPtr_;                            \
                                                                              \
    public: virtual Control& control() const                                  \
    {                                                                         \
        if (controlPtr_.empty())                                              \
        {                                                                     \
            controlPtr_.set(new Control(this->mesh()));                       \
        }                                                                     \
                                                                              \
        return controlPtr_();                                                 \
    }


# define solverManagerRegion_addStorage(TYPE)                                 \
                                                                              \
    solverManagerNestedClass_addVirtual(Storage, storage)


# define solverManagerRegion_end()                                            \
                                                                              \
    solverManagerNestedClass_end()


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class solverManagerRegion Declaration
\*---------------------------------------------------------------------------*/

template <class MESH>
class solverManagerRegion
{
public:

    // Public typedefs

        typedef MESH Mesh;

        typedef solutionControl Control;
        typedef solverManagerStorage<MESH> Storage;


private:

    // Private Data

        //- Reference to args
        const argList& args_;

        //- Reference to time
        Time& time_;

        //- Reference to mesh
        MESH& mesh_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        solverManagerRegion(const solverManagerRegion&);

        //- Disallow default bitwise assignment
        void operator=(const solverManagerRegion&);


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
        virtual void init(const word& init) const = 0;


public:

    // Constructor

        solverManagerRegion
        (
            const argList& args,
            Time& time,
            MESH& mesh
        )
        :
            args_(args), time_(time), mesh_(mesh)
        {};

    // Destructor

        virtual ~solverManagerRegion() {};

    // Member Functions

        //- Control
        virtual Control& control() const = 0;

        //- Storage
        virtual Storage& storage() const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
