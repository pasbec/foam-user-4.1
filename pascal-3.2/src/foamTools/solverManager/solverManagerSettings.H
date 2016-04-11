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
    SM_WRAPPER_NEW(WRAPPER, Settings)


# define SM_REGIONSETTINGS_NEW_INTERFACES()                                   \
                                                                              \
        private: virtual void read() const;


# define SM_REGIONSETTINGS_NEW_INTERFACES_REGION()                            \
                                                                              \
        public: virtual ThisRegion& region() const                            \
        {                                                                     \
            return static_cast<ThisRegion&>                                   \
            (                                                                 \
                manager().regions()[regionLabel]                              \
            );                                                                \
        };                                                                    \
                                                                              \
        public: inline ThisRegion::Control& control() const                   \
        {                                                                     \
            return region().control();                                        \
        };                                                                    \
                                                                              \
        public: inline ThisRegion::Storage& storage() const                   \
        {                                                                     \
            return region().storage();                                        \
        };


# define SM_REGIONSETTINGS_NEW_INTERFACES_GLOBAL()                            \
                                                                              \
        public: inline Manager::Control& control() const                      \
        {                                                                     \
            return manager().control();                                       \
        };                                                                    \
                                                                              \
        public: inline Manager::Storage& storage() const                      \
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
    public solverManagerWrappers::Base<MANAGER, MESH>
{
private:

    // Private typedefs

        typedef solverManagerBase::Settings<MESH> Base;
        typedef solverManagerWrappers::Base<MANAGER, MESH> Wrapper;

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

        Settings(MESH& mesh)
        :
            Base(),
            Wrapper(mesh)
        {};

    // Destructor

        virtual ~Settings() {};

    // Member Functions

        //- Manager
        virtual const MANAGER& manager() const = 0;

        // Basic

            //- Mesh
            virtual MESH& mesh() const
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
private:

    // Private typedefs

        typedef Settings<MANAGER, MESH> Parent;
        typedef solverManagerBase::Settings<MESH> Base;
        typedef solverManagerWrappers::Base<MANAGER, MESH> Wrapper;

        typedef solverManagerBase::Region
        <
            typename MANAGER::ManagerMesh,
            MESH,
            MANAGER::nRegions
        > Region;

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

        RegionSettings(MESH& mesh)
        :
            Parent(mesh)
        {};

    // Destructor

        virtual ~RegionSettings() {};

    // Member Functions

        //- Manager
        virtual const MANAGER& manager() const = 0;

        //- Region
        virtual Region& region() const = 0;
};


/*---------------------------------------------------------------------------*\
            Class solverManagerSettings::GlobalSettings Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER, class MESH>
class GlobalSettings
:
    public Settings<MANAGER, MESH>
{
private:

    // Private typedefs

        typedef Settings<MANAGER, MESH> Parent;
        typedef solverManagerBase::Settings<MESH> Base;
        typedef solverManagerWrappers::Base<MANAGER, MESH> Wrapper;

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

        GlobalSettings(MESH& mesh)
        :
            Parent(mesh)
        {};

    // Destructor

        virtual ~GlobalSettings() {};

    // Member Functions

        //- Manager
        virtual const MANAGER& manager() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solverManagerSettings

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace solverManagerWrappers
{

/*---------------------------------------------------------------------------*\
            Class solverManagerWrappers::RegionSettings Declaration
\*---------------------------------------------------------------------------*/

    template <class MANAGER>
    struct RegionSettings
    {
        typedef typename MANAGER::Mesh Mesh;

        typedef solverManagerSettings::RegionSettings<MANAGER, Mesh> Type;
    };


/*---------------------------------------------------------------------------*\
             Class solverManagerWrappers::GlobalSettings Declaration
\*---------------------------------------------------------------------------*/

    template <class MANAGER>
    struct GlobalSettings
    {
        typedef typename MANAGER::ManagerMesh Mesh;

        typedef solverManagerSettings::GlobalSettings<MANAGER, Mesh> Type;
    };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solverManagerWrappers

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
