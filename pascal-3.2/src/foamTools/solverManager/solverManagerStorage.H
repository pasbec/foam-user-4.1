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
    solverManagerStorage

Description
    ...

\*---------------------------------------------------------------------------*/

#ifndef solverManagerStorage_H
#define solverManagerStorage_H

#include "solverManagerWrapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGIONSTORAGE_NEW_START(WRAPPER)                                  \
                                                                              \
    SM_WRAPPER_NEW(WRAPPER, Storage)


# define SM_REGIONSTORAGE_NEW_INTERFACES()                                    \
                                                                              \
        private: virtual void create(const word& ccase) const;


# define SM_REGIONSTORAGE_NEW_INTERFACES_REGION()                             \
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
        public: inline ThisRegion::Settings& settings() const                 \
        {                                                                     \
            return region().settings();                                       \
        };


# define SM_REGIONSTORAGE_NEW_INTERFACES_GLOBAL()                             \
                                                                              \
        public: inline Manager::Control& control() const                      \
        {                                                                     \
            return manager().control();                                       \
        };                                                                    \
                                                                              \
        public: inline Manager::Settings& settings() const                    \
        {                                                                     \
            return manager().settings();                                      \
        };


# define SM_REGIONSTORAGE_NEW_PUBLIC()                                        \
                                                                              \
        public:


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGIONSTORAGE_NEW()                                               \
                                                                              \
    SM_REGIONSTORAGE_NEW_START(RegionStorage)                                 \
    SM_REGIONSTORAGE_NEW_INTERFACES()                                         \
    SM_REGIONSTORAGE_NEW_INTERFACES_REGION()                                  \
    SM_REGIONSTORAGE_NEW_PUBLIC()


# define SM_GLOBALSTORAGE_NEW()                                               \
                                                                              \
    SM_REGIONSTORAGE_NEW_START(GlobalStorage)                                 \
    SM_REGIONSTORAGE_NEW_INTERFACES()                                         \
    SM_REGIONSTORAGE_NEW_INTERFACES_GLOBAL()                                  \
    SM_REGIONSTORAGE_NEW_PUBLIC()


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGIONSTORAGE_ADD()                                               \
                                                                              \
    SM_WRAPPER_END()                                                          \
    SM_REGION_REGIONSTORAGE_ADD()


# define SM_GLOBALSTORAGE_ADD()                                               \
                                                                              \
    SM_WRAPPER_END()                                                          \
    SM_MANAGER_GLOBALSTORAGE_ADD()


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solverManagerWrappers
{

/*---------------------------------------------------------------------------*\
               Class solverManagerWrappers::Storage Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER, class MESH>
class Storage
:
    public solverManagerBase::Storage<MESH>,
    public solverManagerWrappers::Base<MANAGER, MESH>
{
public:

    // Public typedefs
    
        typedef MESH Mesh;


private:

    // Private typedefs

        typedef solverManagerBase::Storage<MESH> Base;
        typedef solverManagerWrappers::Base<MANAGER, MESH> Wrapper;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        Storage(const Storage&);

        //- Disallow default bitwise assignment
        void operator=(const Storage&);

        // Basic

            //- Settings dictionary
            virtual const dictionary& settingsDict() const = 0;

            //- Storage dictionary
            virtual const dictionary& storageDict() const = 0;

        // Evolution

            //- Create sequence
            virtual void create(const word& ccase) const = 0;


public:

    // Constructor

        Storage(MESH& mesh)
        :
            Base(),
            Wrapper(mesh)
        {};

    // Destructor

        virtual ~Storage() {};

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
                return storageDict();
            };
};


/*---------------------------------------------------------------------------*\
              Class solverManagerStorage::RegionStorage Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER>
class RegionStorage
:
    public Storage<MANAGER, typename MANAGER::Mesh>
{
public:

    // Public typedefs

        typedef typename MANAGER::Mesh Mesh;


private:

    // Private typedefs

        typedef Storage<MANAGER, Mesh> Parent;
        typedef solverManagerBase::Storage<Mesh> Base;
        typedef solverManagerWrappers::Base<MANAGER, Mesh> Wrapper;

        typedef solverManagerBase::Region_
        <
            typename MANAGER::ManagerMesh,
            Mesh,
            MANAGER::nRegions
        > Region_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        RegionStorage(const RegionStorage&);

        //- Disallow default bitwise assignment
        void operator=(const RegionStorage&);

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

            //- Create sequence
            virtual void create(const word& ccase) const = 0;


public:

    // Constructor

        RegionStorage(Mesh& mesh)
        :
            Parent(mesh)
        {};

    // Destructor

        virtual ~RegionStorage() {};

    // Member Functions

        //- Manager
        virtual const MANAGER& manager() const = 0;

        //- Region
        virtual Region_& region() const = 0;
};


/*---------------------------------------------------------------------------*\
             Class solverManagerStorage::GlobalStorage Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER>
class GlobalStorage
:
    public Storage<MANAGER, typename MANAGER::ManagerMesh>
{
public:

    // Public typedefs

        typedef typename MANAGER::ManagerMesh Mesh;


private:

    // Private typedefs

        typedef Storage<MANAGER, Mesh> Parent;
        typedef solverManagerBase::Storage<Mesh> Base;
        typedef solverManagerWrappers::Base<MANAGER, Mesh> Wrapper;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        GlobalStorage(const GlobalStorage&);

        //- Disallow default bitwise assignment
        void operator=(const GlobalStorage&);

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

            //- Create sequence
            virtual void create(const word& ccase) const = 0;


public:

    // Constructor

        GlobalStorage(Mesh& mesh)
        :
            Parent(mesh)
        {};

    // Destructor

        virtual ~GlobalStorage() {};

    // Member Functions

        //- Manager
        virtual const MANAGER& manager() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solverManagerWrappers

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
