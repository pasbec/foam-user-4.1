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

#include "solverManager.H"
#include "solverManagerWrapper.H"
#include "solverManagerRegion.H"
#include "solverManagerStorage.H"

// TODO: Get the correct mesh in the Wrapper! Therefor I have to change:
//       public: typedef typename MANAGER::ManagerMesh Mesh;
//       But ONLY FOR REGION, SETTINGS, and STORAGE of Regions!!!

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define smRegions_new(MANAGER, NAME)                                         \
                                                                              \
    smWrapper_new(Regions, solverManagerRegions, MANAGER, NAME)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define smRegions_end()                                                      \
                                                                              \
    smWrapper_end()


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define smRegions_addRegion_manager(TYPE, NAME)                              \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: virtual TYPE& NAME() const                                    \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                smWrapper_setItem_P2(TYPE, NAME, this->manager(), this->mesh())\
                                                                              \
                regions().insert                                              \
                (                                                             \
                    NAME##Ptr_->mesh().name(),                                \
                    NAME##Ptr_.operator->()                                   \
                );                                                            \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };

# define smRegions_addDefaultRegion_manager(TYPE)                             \
                                                                              \
        smRegions_addRegion_manager(TYPE, defaultRegion)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class solverManagerRegions Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER, class MESH>
class solverManagerRegions
:
    public solverManagerWrapper<MANAGER, MESH>
{
public:

    // Public typedefs

        typedef MANAGER Manager;
        typedef MESH Mesh;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        solverManagerRegions(const solverManagerRegions&);

        //- Disallow default bitwise assignment
        void operator=(const solverManagerRegions&);


public:

    // Constructor

        solverManagerRegions(const Manager& manager, Mesh& mesh)
        :
            solverManagerWrapper<Manager, Mesh>
            (
                manager,
                mesh
            )
        {};

    // Destructor

        virtual ~solverManagerRegions() {};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
