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
    solverManagerSettings

Description
    ...

\*---------------------------------------------------------------------------*/

#ifndef solverManagerSettings_H
#define solverManagerSettings_H

#include "solverManagerRegion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_SETTINGS_NEW(MANAGER)                                             \
                                                                              \
    SM_WRAPPER_NEW(Settings, MANAGER, Settings)                               \
                                                                              \
    private: virtual void read() const;                                       \
                                                                              \
    public:


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_SETTINGS_MEMBER_ADD(TYPE, NAME)                                   \
                                                                              \
        public: mutable TYPE NAME;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_SETTINGS_ADD()                                                    \
                                                                              \
    SM_WRAPPER_END()                                                          \
    SM_REGION_SETTINGS_ADD()


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class solverManagerSettings Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER>
class solverManagerSettings
:
    public solverManagerRegion<MANAGER>
{
public:

    // Public typedefs

        typedef MANAGER Manager;
        typedef typename MANAGER::Settings Base;

        typedef solverManagerRegion<MANAGER> Wrapper;

        typedef typename Wrapper::ManagerMesh ManagerMesh;
        typedef typename Wrapper::Mesh Mesh;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        solverManagerSettings(const solverManagerSettings&);

        //- Disallow default bitwise assignment
        void operator=(const solverManagerSettings&);


public:

    // Constructor

        solverManagerSettings(const Manager& manager, Mesh& mesh)
        :
            Wrapper(manager, mesh)
        {};

    // Destructor

        virtual ~solverManagerSettings() {};

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
