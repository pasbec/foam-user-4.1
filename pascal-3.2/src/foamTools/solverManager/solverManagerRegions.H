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

#include "solverManagerWrapper.H"

// TODO: Fix the invalid pointer bug during destruction!

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGIONS_NEW(MANAGER)                                              \
                                                                              \
    SM_WRAPPER_NEW(Regions, MANAGER, Regions)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGIONS_REGION_ADD(TYPE, NAME, PMESH)                             \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: virtual TYPE& NAME() const                                    \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                SM_WRAPPER_MEMBER_SETINIT_P2(                                 \
                    TYPE, NAME, this->manager(), PMESH)                       \
                                                                              \
                this->Base::regions().insert                                  \
                (                                                             \
                    NAME##Ptr_->mesh().name(),                                \
                    NAME##Ptr_.operator->()                                   \
                );                                                            \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };

# define SM_REGIONS_DEFAULTREGION_ADD(PMESH)                                  \
                                                                              \
        SM_REGIONS_REGION_ADD(DefaultRegion, defaultRegion, PMESH)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGIONS_ADD()                                                     \
                                                                              \
    SM_WRAPPER_END()                                                          \
    SM_ADD_REGIONS()


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class solverManagerRegions Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER>
class solverManagerRegions
:
    public solverManagerWrapper<MANAGER, typename MANAGER::ManagerMesh>
{
public:

    // Public typedefs

        typedef MANAGER Manager;
        typedef typename MANAGER::Regions Base;

        typedef typename MANAGER::ManagerMesh ManagerMesh;
        typedef typename MANAGER::ManagerMesh Mesh;

        typedef solverManagerWrapper<Manager, Mesh> Wrapper;


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
            Wrapper(manager, mesh)
        {};

    // Destructor

        virtual ~solverManagerRegions() {};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
