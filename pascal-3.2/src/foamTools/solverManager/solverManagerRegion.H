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

# define SM_REGION_CONTROL_ADD(TYPE)                                          \
                                                                              \
        public: typedef TYPE Control;                                         \
                                                                              \
        SM_WRAPPER_MEMBER_MOD_NOINIT_ADD_P1(                                  \
            Control, control, virtual, this->mesh())


# define SM_REGION_SETTINGS_ADD()                                             \
                                                                              \
        SM_WRAPPER_MEMBER_VIRTUAL_ADD_P2(                                     \
            Settings, settings, this->manager(), this->mesh())


# define SM_REGION_STORAGE_ADD()                                              \
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

template <class MANAGER>
class solverManagerRegion
:
    public solverManagerWrapper<MANAGER, typename MANAGER::Mesh>
{
public:

    // Public typedefs

        typedef MANAGER Manager;
        typedef typename MANAGER::Region Base;

        typedef typename MANAGER::ManagerMesh ManagerMesh;
        typedef typename MANAGER::Mesh Mesh;

        typedef solverManagerWrapper<Manager, Mesh> Wrapper;


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
            Wrapper(manager, mesh)
        {};

    // Destructor

        virtual ~solverManagerRegion() {};

    // Member Functions

        // Basic

            //- Reference to settings dictionary
            inline const dictionary& settingsDict() const
            {
                return this->Wrapper::settingsDict().subDict
                (
                    this->Wrapper::mesh().name()
                );
            };

            //- Reference to storage dictionary
            inline const dictionary& storageDict() const
            {
                return this->Wrapper::storageDict().subDict
                (
                    this->Wrapper::mesh().name()
                );
            };
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
