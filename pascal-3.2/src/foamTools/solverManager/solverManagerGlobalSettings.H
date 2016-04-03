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
    solverManagerGlobalSettings

Description
    ...

\*---------------------------------------------------------------------------*/

#ifndef solverManagerGlobalSettings_H
#define solverManagerGlobalSettings_H

#include "solverManagerWrapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_GLOBALSETTINGS_NEW(MANAGER)                                       \
                                                                              \
    SM_WRAPPER_NEW(GlobalSettings, MANAGER, GlobalSettings)                   \
                                                                              \
    private: virtual void read() const;                                       \
                                                                              \
    public:


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_GLOBALSETTINGS_MEMBER_ADD(TYPE, NAME)                             \
                                                                              \
        public: mutable TYPE NAME;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_GLOBALSETTINGS_ADD()                                              \
                                                                              \
    SM_WRAPPER_END()                                                          \
    SM_ADD_GLOBALSETTINGS()


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class solverManagerGlobalSettings Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER>
class solverManagerGlobalSettings
:
    public solverManagerWrapper<MANAGER, typename MANAGER::ManagerMesh>
{
public:

    // Public typedefs

        typedef MANAGER Manager;
        typedef typename MANAGER::GlobalSettings Base;

        typedef typename MANAGER::ManagerMesh ManagerMesh;
        typedef typename MANAGER::ManagerMesh Mesh;

        typedef solverManagerWrapper<Manager, Mesh> Wrapper;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        solverManagerGlobalSettings(const solverManagerGlobalSettings&);

        //- Disallow default bitwise assignment
        void operator=(const solverManagerGlobalSettings&);


public:

    // Constructor

        solverManagerGlobalSettings(const Manager& manager, Mesh& mesh)
        :
            solverManagerWrapper<Manager, Mesh>
            (
                manager,
                mesh
            )
        {};

    // Destructor

        virtual ~solverManagerGlobalSettings() {};

    // Member Functions

        // Basic

            //- Reference to dictionary
            inline const dictionary& dict() const
            {
                return this->Wrapper::settingsDict();
            };
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
