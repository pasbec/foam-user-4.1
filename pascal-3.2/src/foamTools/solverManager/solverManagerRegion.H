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

#include "solverManagerWrapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGION_NEW(MANAGER, NAME)                                         \
                                                                              \
    SM_WRAPPER_NEW(Region, MANAGER, NAME)

# define SM_DEFAULTREGION_NEW(MANAGER)                                        \
                                                                              \
    SM_REGION_NEW(MANAGER, DefaultRegion)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGION_REGIONSETTINGS_ADD()                                       \
                                                                              \
        SM_WRAPPER_MEMBER_VIRTUAL_ADD_P2(                                     \
            Settings, settings, this->manager(), this->mesh())


# define SM_REGION_REGIONSTORAGE_ADD()                                        \
                                                                              \
        SM_WRAPPER_MEMBER_VIRTUAL_ADD_P2(                                     \
            Storage, storage, this->manager(), this->mesh())


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGION_END()                                                      \
                                                                              \
    SM_WRAPPER_END()


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class solverManagerRegion Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER, class MESH>
class solverManagerRegion
:
    public solverManagerBase::Region<MESH>,
    public solverManagerWrapper<MANAGER, MESH>
{
public:

    // Public typedefs

        typedef MANAGER Manager;

        typedef typename Manager::ManagerMesh ManagerMesh;
        typedef MESH Mesh;

        typedef solverManagerBase::Region<Mesh> Base;
        typedef solverManagerWrapper<Manager, Mesh> Wrapper;

        typedef typename Manager::Control Control;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        solverManagerRegion(const solverManagerRegion&);

        //- Disallow default bitwise assignment
        void operator=(const solverManagerRegion&);


public:

    // Constructor

        solverManagerRegion(const Manager& manager, Mesh& mesh)
        :
            Base(),
            Wrapper(manager, mesh)
        {};

    // Destructor

        virtual ~solverManagerRegion() {};

    // Member Functions

        //- Control
        virtual Control& control() const
        {
            return this->manager().control();
        };

        // Basic

            //- Mesh
            inline Mesh& mesh() const
            {
                return this->Wrapper::wrappedMesh();
            }

        // Dictionaries

            //- Settings dictionary
            inline const dictionary& settingsDict() const
            {
                return this->Wrapper::settingsDict().subDict
                (
                    this->Wrapper::mesh().name()
                );
            };

            //- Storage dictionary
            inline const dictionary& storageDict() const
            {
                return this->Wrapper::storageDict().subDict
                (
                    this->Wrapper::mesh().name()
                );
            };
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solverManagerWrappers
{
    template <class MANAGER>
    struct Region
    {
        typedef MANAGER Manager;

        typedef typename MANAGER::ManagerMesh ManagerMesh;
        typedef typename MANAGER::Mesh Mesh;

        typedef solverManagerRegion<Manager, Mesh> Type;
    };
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
