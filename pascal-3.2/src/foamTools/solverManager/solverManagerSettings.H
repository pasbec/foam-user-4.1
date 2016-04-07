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

#include "solverManagerWrapper.H"
#include "solverManagerRegion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGIONSETTINGS_NEW(MANAGER)                                       \
                                                                              \
    SM_WRAPPER_NEW(RegionSettings, MANAGER, Settings)                         \
                                                                              \
    private: virtual void read() const;                                       \
                                                                              \
    public:


# define SM_GLOBALSETTINGS_NEW(MANAGER)                                       \
                                                                              \
    SM_WRAPPER_NEW(GlobalSettings, MANAGER, Settings)                         \
                                                                              \
    private: virtual void read() const;                                       \
                                                                              \
    public:


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_SETTINGS_MEMBER_ADD(TYPE, NAME)                                   \
                                                                              \
        public: mutable TYPE NAME;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGIONSETTINGS_ADD()                                              \
                                                                              \
    SM_WRAPPER_END()                                                          \
    SM_REGION_REGIONSETTINGS_ADD()


# define SM_GLOBALSETTINGS_ADD()                                              \
                                                                              \
    SM_WRAPPER_END()                                                          \
    SM_MANAGER_GLOBALSETTINGS_ADD()


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solverManagerSettings
{

/*---------------------------------------------------------------------------*\
                       Class RegionSettings Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER, class MESH>
class RegionSettings
:
    public solverManagerBase::Settings<MESH>,
    public solverManagerWrapper<MANAGER, MESH>
{
public:

    // Public typedefs

        typedef MANAGER Manager;

        typedef typename Manager::ManagerMesh ManagerMesh;
        typedef MESH Mesh;

        typedef solverManagerBase::Settings<Mesh> Base;
        typedef solverManagerWrapper<Manager, Mesh> Wrapper;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        RegionSettings(const RegionSettings&);

        //- Disallow default bitwise assignment
        void operator=(const RegionSettings&);

        // Basic

            //- Reference to settings dictionary
            inline const dictionary& settingsDict() const
            {
                return this->Wrapper::settingsDict().subDict
                (
                    mesh().name()
                );
            };

            //- Reference to storage dictionary
            inline const dictionary& storageDict() const
            {
                return this->Wrapper::storageDict().subDict
                (
                    mesh().name()
                );
            };

            
public:

    // Constructor

        RegionSettings(const Manager& manager, Mesh& mesh)
        :
            Base(),
            Wrapper(manager, mesh)
        {};

    // Destructor

        virtual ~RegionSettings() {};

    // Member Functions

        // Basic

            //- Mesh
            inline Mesh& mesh() const
            {
                return this->Wrapper::wrappedMesh();
            }

        // Dictionaries

            //- Dictionary
            inline const dictionary& dict() const
            {
                return settingsDict();
            };
};


/*---------------------------------------------------------------------------*\
                       Class GlobalSettings Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER, class MESH>
class GlobalSettings
:
    public solverManagerBase::Settings<MESH>,
    public solverManagerWrapper<MANAGER, MESH>
{
public:

    // Public typedefs

        typedef MANAGER Manager;

        typedef typename Manager::ManagerMesh ManagerMesh;
        typedef MESH Mesh;

        typedef solverManagerBase::Settings<Mesh> Base;
        typedef solverManagerWrapper<Manager, Mesh> Wrapper;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        GlobalSettings(const GlobalSettings&);

        //- Disallow default bitwise assignment
        void operator=(const GlobalSettings&);

        // Basic

            //- Reference to settings dictionary
            inline const dictionary& settingsDict() const
            {
                return this->Wrapper::settingsDict();
            };

            //- Reference to storage dictionary
            inline const dictionary& storageDict() const
            {
                return this->Wrapper::storageDict();
            };

            
public:

    // Constructor

        GlobalSettings(const Manager& manager, Mesh& mesh)
        :
            Base(),
            Wrapper(manager, mesh)
        {};

    // Destructor

        virtual ~GlobalSettings() {};

    // Member Functions

        // Basic

            //- Mesh
            inline Mesh& mesh() const
            {
                return this->Wrapper::wrappedMesh();
            }

        // Dictionaries

            //- Dictionary
            inline const dictionary& dict() const
            {
                return settingsDict();
            };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solverManagerSettings

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solverManagerWrappers
{

    template <class MANAGER>
    struct RegionSettings
    {
        typedef MANAGER Manager;

        typedef typename MANAGER::ManagerMesh ManagerMesh;
        typedef typename MANAGER::Mesh Mesh;

        typedef solverManagerSettings::RegionSettings<Manager, Mesh> Type;
    };
    
    template <class MANAGER>
    struct GlobalSettings
    {
        typedef MANAGER Manager;

        typedef typename MANAGER::ManagerMesh ManagerMesh;
        typedef typename MANAGER::ManagerMesh Mesh;

        typedef solverManagerSettings::GlobalSettings<Manager, Mesh> Type;
    };
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
