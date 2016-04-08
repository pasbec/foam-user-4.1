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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGIONS_NEW()                                                     \
                                                                              \
    SM_WRAPPER_NEW(Regions, Regions)                                          \
                                                                              \
        public: typedef Regions RegionsType;                                  \


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGIONS_REGION_ADD(TYPE, NAME, MESH)                              \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: virtual TYPE& NAME() const                                    \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                SM_WRAPPER_MEMBER_SETINIT_P2(                                 \
                    TYPE, NAME, this->manager(), MESH)                        \
                                                                              \
                this->Base::regionHashTable().insert                          \
                (                                                             \
                    NAME##Ptr_->mesh().name(),                                \
                    NAME##Ptr_.operator->()                                   \
                );                                                            \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };


# define SM_REGIONS_DEFAULTREGION_ADD(MESH)                                   \
                                                                              \
        SM_REGIONS_REGION_ADD(DefaultRegion, defaultRegion, MESH)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGIONS_ADD()                                                     \
                                                                              \
    SM_WRAPPER_END()                                                          \
    SM_MANAGER_REGIONS_ADD()


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class solverManagerRegions Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER>
class solverManagerRegions
:
    public solverManagerBase::Regions
    <
        typename MANAGER::ManagerMesh, typename MANAGER::Mesh
    >,
    public solverManagerWrapper<MANAGER, typename MANAGER::ManagerMesh>
{
public:

    // Public typedefs

        typedef MANAGER Manager;

        typedef typename Manager::ManagerMesh ManagerMesh;
        typedef typename Manager::ManagerMesh Mesh;
        typedef typename Manager::Mesh SubMesh;

        typedef solverManagerBase::Regions
        <
            typename Manager::ManagerMesh, SubMesh
        > Base;
        typedef solverManagerWrapper<Manager, Mesh> Wrapper;

        typedef typename solverManagerBase::Region<SubMesh> Region;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        solverManagerRegions(const solverManagerRegions&);

        //- Disallow default bitwise assignment
        void operator=(const solverManagerRegions&);


public:

    // Constructor

        solverManagerRegions(Mesh& mesh)
        :
            Base(),
            Wrapper(mesh)
        {};

    // Destructor

        virtual ~solverManagerRegions() {};

    // Member Functions

        //- Manager
        virtual const Manager& manager() const = 0;

        // Basic

            //- Mesh
            virtual Mesh& mesh() const
            {
                return this->Wrapper::wrappedMesh();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solverManagerWrappers
{
    template <class MANAGER>
    struct Regions
    {
        typedef MANAGER Manager;

        typedef typename MANAGER::ManagerMesh ManagerMesh;
        typedef typename MANAGER::ManagerMesh Mesh;

        typedef solverManagerRegions<Manager> Type;
    };
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
