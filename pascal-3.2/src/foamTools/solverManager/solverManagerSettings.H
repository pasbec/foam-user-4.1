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

# define SM_REGIONSETTINGS_NEW_START(WRAPPER)                                 \
                                                                              \
    SM_WRAPPER_NEW(WRAPPER, Settings)                                         \
                                                                              \
        public: typedef Settings SettingsType;


# define SM_REGIONSETTINGS_NEW_INTERFACES()                                   \
                                                                              \
        private: virtual void read() const;


# define SM_REGIONSETTINGS_NEW_INTERFACES_REGION()                            \
                                                                              \
        public: virtual RegionType& region() const                            \
        {                                                                     \
            return static_cast<RegionType&>                                   \
            (                                                                 \
                manager().regions()[this->mesh().name()]                      \
            );                                                                \
        };                                                                    \
                                                                              \
        public: virtual RegionType::Control& control() const                  \
        {                                                                     \
            return region().control();                                        \
        };                                                                    \
                                                                              \
        public: virtual RegionType::Storage& storage() const                  \
        {                                                                     \
            return region().storage();                                        \
        };


# define SM_REGIONSETTINGS_NEW_INTERFACES_GLOBAL()                            \
                                                                              \
        public: virtual ManagerType::Control& control() const                 \
        {                                                                     \
            return manager().control();                                       \
        };                                                                    \
                                                                              \
        public: virtual ManagerType::Storage& storage() const                 \
        {                                                                     \
            return manager().storage();                                       \
        };



# define SM_REGIONSETTINGS_NEW_PUBLIC()                                       \
                                                                              \
        public:


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGIONSETTINGS_NEW()                                              \
                                                                              \
    SM_REGIONSETTINGS_NEW_START(RegionSettings)                               \
    SM_REGIONSETTINGS_NEW_INTERFACES()                                        \
    SM_REGIONSETTINGS_NEW_INTERFACES_REGION()                                 \
    SM_REGIONSETTINGS_NEW_PUBLIC()


# define SM_GLOBALSETTINGS_NEW()                                              \
                                                                              \
    SM_REGIONSETTINGS_NEW_START(GlobalSettings)                               \
    SM_REGIONSETTINGS_NEW_INTERFACES()                                        \
    SM_REGIONSETTINGS_NEW_INTERFACES_GLOBAL()                                 \
    SM_REGIONSETTINGS_NEW_PUBLIC()


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
                Class solverManagerSettings::Settings Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER, class MESH>
class Settings
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

        typedef typename solverManagerBase::Control Control;
        typedef typename solverManagerBase::Storage<Mesh> Storage;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        Settings(const Settings&);

        //- Disallow default bitwise assignment
        void operator=(const Settings&);

        // Basic

            //- Settings dictionary
            virtual const dictionary& settingsDict() const = 0;;

            //- Storage dictionary
            virtual const dictionary& storageDict() const = 0;

        // Evolution

            //- Read sequence
            virtual void read() const = 0;


public:

    // Constructor

        Settings(Mesh& mesh)
        :
            Base(),
            Wrapper(mesh)
        {};

    // Destructor

        virtual ~Settings() {};

    // Member Functions

        //- Manager
        virtual const Manager& manager() const = 0;

        //- Control
        virtual Control& control() const = 0;

        //- Storage
        virtual Storage& storage() const = 0;

        // Basic

            //- Mesh
            virtual Mesh& mesh() const
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
             Class solverManagerSettings::RegionSettings Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER, class MESH>
class RegionSettings
:
    public Settings<MANAGER, MESH>
{
public:

    // Public typedefs

        typedef MANAGER Manager;

        typedef typename Manager::ManagerMesh ManagerMesh;
        typedef MESH Mesh;

        typedef Settings<MANAGER, MESH> Parent;
        typedef typename Parent::Base Base;
        typedef typename Parent::Wrapper Wrapper;

        typedef typename solverManagerBase::Region<Mesh> Region;
        typedef typename Parent::Control Control;
        typedef typename Parent::Storage Storage;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        RegionSettings(const RegionSettings&);

        //- Disallow default bitwise assignment
        void operator=(const RegionSettings&);

        // Basic

            //- Settings dictionary
            virtual const dictionary& settingsDict() const
            {
                return this->Wrapper::settingsDict().subDict
                (
                    this->Parent::mesh().name()
                );
            };

            //- Storage dictionary
            virtual const dictionary& storageDict() const
            {
                return this->Wrapper::storageDict().subDict
                (
                    this->Parent::mesh().name()
                );
            };

        // Evolution

            //- Read sequence
            virtual void read() const = 0;


public:

    // Constructor

        RegionSettings(Mesh& mesh)
        :
            Parent(mesh)
        {};

    // Destructor

        virtual ~RegionSettings() {};

    // Member Functions

        //- Manager
        virtual const Manager& manager() const = 0;

        //- Region
        virtual Region& region() const = 0;

        //- Control
        virtual Control& control() const = 0;

        //- Storage
        virtual Storage& storage() const = 0;
};


/*---------------------------------------------------------------------------*\
            Class solverManagerSettings::GlobalSettings Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER, class MESH>
class GlobalSettings
:
    public Settings<MANAGER, MESH>
{
public:

    // Public typedefs

        typedef MANAGER Manager;

        typedef typename Manager::ManagerMesh ManagerMesh;
        typedef MESH Mesh;

        typedef Settings<MANAGER, MESH> Parent;
        typedef typename Parent::Base Base;
        typedef typename Parent::Wrapper Wrapper;
 
        typedef typename Parent::Control Control;
        typedef typename Parent::Storage Storage;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        GlobalSettings(const GlobalSettings&);

        //- Disallow default bitwise assignment
        void operator=(const GlobalSettings&);

        // Basic

            //- Settings dictionary
            virtual const dictionary& settingsDict() const
            {
                return this->Wrapper::settingsDict();
            };

            //- Storage dictionary
            virtual const dictionary& storageDict() const
            {
                return this->Wrapper::storageDict();
            };

        // Evolution

            //- Read sequence
            virtual void read() const = 0;


public:

    // Constructor

        GlobalSettings(Mesh& mesh)
        :
            Parent(mesh)
        {};

    // Destructor

        virtual ~GlobalSettings() {};

    // Member Functions

        //- Manager
        virtual const Manager& manager() const = 0;

        //- Control
        virtual Control& control() const = 0;

        //- Storage
        virtual Storage& storage() const = 0;
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
