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
#include "solverManagerRegion.H"

// TODO: Remove old storage

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGIONSTORAGE_NEW_START(WRAPPER)                                  \
                                                                              \
    SM_WRAPPER_NEW(WRAPPER, Storage)                                          \
                                                                              \
        public: typedef Storage StorageType;


# define SM_REGIONSTORAGE_NEW_INTERFACES()                                    \
                                                                              \
        private: virtual void create() const;


# define SM_REGIONSTORAGE_NEW_INTERFACES_REGION()                             \
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
        public: virtual RegionType::Settings& settings() const                \
        {                                                                     \
            return region().settings();                                       \
        };


# define SM_REGIONSTORAGE_NEW_INTERFACES_GLOBAL()                             \
                                                                              \
        public: virtual ManagerType::Control& control() const                 \
        {                                                                     \
            return manager().control();                                       \
        };                                                                    \
                                                                              \
        public: virtual ManagerType::Settings& settings() const               \
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

# define SM_STORAGE_MEMBER_ADD(TYPE, NAME)                                    \
                                                                              \
    private: mutable autoPtr<TYPE> NAME##Ptr_;                                \
                                                                              \
    private: mutable Switch NAME##Optional_;                                  \
                                                                              \
    private: mutable Switch NAME##Locked_;                                    \
                                                                              \
    private: mutable Switch NAME##Disabled_;                                  \
                                                                              \
    private: inline void mandatory_##NAME() const                             \
    {                                                                         \
        NAME##Optional_ = false;                                              \
    };                                                                        \
                                                                              \
    private: inline void optional_##NAME() const                              \
    {                                                                         \
        NAME##Optional_ = true;                                               \
    };                                                                        \
                                                                              \
    private: void enable_##NAME() const                                       \
    {                                                                         \
        if (NAME##Disabled_)                                                  \
        {                                                                     \
            NAME##Disabled_ = false;                                          \
        }                                                                     \
        else                                                                  \
        {                                                                     \
            FatalErrorIn                                                      \
            (this->manager().type() + "Storage::enable_"+#NAME+"()")          \
                << "The storage item "                                        \
                    << #NAME                                                  \
                    << " is already disabled!"                                \
                    << abort(FatalError);                                     \
        }                                                                     \
    };                                                                        \
                                                                              \
    private: void disable_##NAME() const                                      \
    {                                                                         \
        if (!NAME##Disabled_)                                                 \
        {                                                                     \
            if (NAME##Ptr_.valid())                                           \
            {                                                                 \
                FatalErrorIn                                                  \
                (this->manager().type() + "Storage::disable_"+#NAME+"()")     \
                    << "The storage item "                                    \
                        << #NAME                                              \
                        << " is still allocated and cannot be disabled."      \
                        << " Clear it first!"                                 \
                        << abort(FatalError);                                 \
                                                                              \
            }                                                                 \
            else                                                              \
            {                                                                 \
                NAME##Disabled_ = true;                                       \
            }                                                                 \
        }                                                                     \
    };                                                                        \
                                                                              \
    private: void create_##NAME(const word& init = "default") const           \
    {                                                                         \
        if (!NAME##Disabled_)                                                 \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                if (!NAME##Locked_)                                           \
                {                                                             \
                    NAME##Locked_ = true;                                     \
                                                                              \
                    Info << "Create " << #NAME                                \
                        << " for region " << mesh().name()                    \
                        << endl;                                              \
                                                                              \
                    init_##NAME(init);                                        \
                                                                              \
                    NAME##Locked_ = false;                                    \
                }                                                             \
                else                                                          \
                {                                                             \
                    FatalErrorIn                                              \
                    (this->manager().type() + "Storage::create_"+#NAME+"()")  \
                        << "Pointer for "                                     \
                            << #NAME                                          \
                            << " is part of a circular reference."            \
                            << " Check all occurances in storage!"            \
                            << abort(FatalError);                             \
                }                                                             \
            }                                                                 \
                                                                              \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                FatalErrorIn                                                  \
                (this->manager().type() + "Storage::create_"+#NAME+"()")      \
                    << "Pointer for "                                         \
                        << #NAME                                              \
                        << " is still unallocated after creation."            \
                        << " Make sure the data member " << #NAME << "Ptr_ "  \
                        << " is set correctly in its init-function: "         \
                        << "init_" << #NAME << "(const word& init)!"          \
                        << abort(FatalError);                                 \
            }                                                                 \
        }                                                                     \
        else                                                                  \
        {                                                                     \
            FatalErrorIn                                                      \
            (this->manager().type() + "Storage::create_"+#NAME+"()")          \
                << "The storage item "                                        \
                    << #NAME                                                  \
                    << " is globally disabled and cannot be used."            \
                    << " Check your storage settings!"                        \
                    << abort(FatalError);                                     \
                                                                              \
        }                                                                     \
    };                                                                        \
                                                                              \
    private: bool dict_##NAME                                                 \
    (                                                                         \
        dictionary& dict,                                                     \
        const word& name = #NAME                                              \
    ) const                                                                   \
    {                                                                         \
        bool isValid = itemDict(name, dict);                                  \
                                                                              \
        if (!NAME##Optional_)                                                 \
        {                                                                     \
            return !isValid;                                                  \
        }                                                                     \
        else                                                                  \
        {                                                                     \
            return isValid;                                                   \
        }                                                                     \
    };                                                                        \
                                                                              \
    private: bool switch_##NAME                                               \
    (                                                                         \
        dictionary& dict,                                                     \
        const word& name = #NAME                                              \
    ) const                                                                   \
    {                                                                         \
        disable_##NAME();                                                     \
                                                                              \
        if (dict_##NAME(dict, name))                                          \
        {                                                                     \
            enable_##NAME();                                                  \
                                                                              \
            if (NAME##Optional_)                                              \
            {                                                                 \
                Info << "Enable " << #NAME                                    \
                    << " for region " << mesh().name()                        \
                    << endl;                                                  \
            }                                                                 \
                                                                              \
            return true;                                                      \
        }                                                                     \
                                                                              \
        if (!NAME##Optional_)                                                 \
        {                                                                     \
            FatalErrorIn                                                      \
            (this->manager().type() + "Storage::switch_"+#NAME+"()")          \
                        << "Storage element "                                 \
                            << #NAME                                          \
                            << " is mandatory and cannot be disabled."        \
                            << " Check settings in storage and properties!"   \
                            << abort(FatalError);                             \
        }                                                                     \
        else                                                                  \
        {                                                                     \
            Info << "Disable " << #NAME                                       \
                << " for region " << mesh().name()                            \
                << endl;                                                      \
        }                                                                     \
                                                                              \
        return false;                                                         \
    };                                                                        \
                                                                              \
    private: bool make_##NAME(const word& name = #NAME) const                 \
    {                                                                         \
        word init = "";                                                       \
                                                                              \
        dictionary dict;                                                      \
        if (switch_##NAME(dict, name))                                        \
        {                                                                     \
            init = dict.lookupOrDefault("init", word("default"));             \
                                                                              \
            create_##NAME(init);                                              \
                                                                              \
            return true;                                                      \
        }                                                                     \
                                                                              \
        return false;                                                         \
    };                                                                        \
                                                                              \
    private: bool makeMandatory_##NAME(const word& name = #NAME) const        \
    {                                                                         \
         mandatory_##NAME();                                                  \
                                                                              \
         return make_##NAME(name);                                            \
    };                                                                        \
                                                                              \
    private: bool makeOptional_##NAME(const word& name = #NAME) const         \
    {                                                                         \
         optional_##NAME();                                                   \
                                                                              \
         return make_##NAME(name);                                            \
    };                                                                        \
                                                                              \
    private: void init_##NAME(const word& init = "default") const;            \
                                                                              \
    public: inline bool is_##NAME() const                                     \
    {                                                                         \
        return !NAME##Disabled_;                                              \
    };                                                                        \
                                                                              \
    public: inline bool isNot_##NAME() const                                  \
    {                                                                         \
        return NAME##Disabled_;                                               \
    };                                                                        \
                                                                              \
    public: inline bool isEmpty_##NAME() const                                \
    {                                                                         \
        return NAME##Ptr_.empty();                                            \
    };                                                                        \
                                                                              \
    public: inline bool isValid_##NAME() const                                \
    {                                                                         \
        return NAME##Ptr_.valid();                                            \
    };                                                                        \
                                                                              \
    public: TYPE& NAME(const word& init = "default") const                    \
    {                                                                         \
        create_##NAME(init);                                                  \
                                                                              \
        return NAME##Ptr_();                                                  \
    };                                                                        \
                                                                              \
    public: TYPE* NAME##Ptr(const word& init = "default") const               \
    {                                                                         \
        create_##NAME(init);                                                  \
                                                                              \
        return NAME##Ptr_.operator->();                                       \
    };


# define SM_GLOBALSTORAGE_MEMBER_ADD(TYPE, NAME)                              \
                                                                              \
        SM_STORAGE_MEMBER_ADD(TYPE, NAME)


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

namespace solverManagerStorage
{

/*---------------------------------------------------------------------------*\
              Class solverManagerStorage::Storage Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER, class MESH>
class Storage
:
    public solverManagerBase::Storage<MESH>,
    public solverManagerWrapper<MANAGER, MESH>
{
public:

    // Public typedefs

        typedef MANAGER Manager;

        typedef typename Manager::ManagerMesh ManagerMesh;
        typedef MESH Mesh;

        typedef solverManagerBase::Storage<Mesh> Base;
        typedef solverManagerWrapper<Manager, Mesh> Wrapper;

        typedef typename solverManagerBase::Control Control;
        typedef typename solverManagerBase::Settings<Mesh> Settings;


private:

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
            virtual void create() const = 0;


public:

    // Constructor

        Storage(Mesh& mesh)
        :
            Base(),
            Wrapper(mesh)
        {};

    // Destructor

        virtual ~Storage() {};

    // Member Functions

        //- Manager
        virtual const Manager& manager() const = 0;

        //- Control
        virtual Control& control() const = 0;

        //- Settings
        virtual Settings& settings() const = 0;

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
                return storageDict();
            };


// TODO: Remove after old storage is removed
            //- Read dictionary of named item in storage (if present)
            bool itemDict
            (
                const word& name,
                dictionary& itemDict
            ) const
            {
                if
                (
                    dict().found(name)
                 && dict().isDict(name)
                )
                {
                    itemDict = dict().subDict(name);

                    Switch enabled;

                    if (itemDict.readIfPresent("enabled", enabled))
                    {
                        return enabled;
                    }
                }

                return false;
            };
};


/*---------------------------------------------------------------------------*\
              Class solverManagerStorage::RegionStorage Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER, class MESH>
class RegionStorage
:
    public Storage<MANAGER, MESH>
{
public:

    // Public typedefs

        typedef MANAGER Manager;

        typedef typename Manager::ManagerMesh ManagerMesh;
        typedef MESH Mesh;

        typedef Storage<MANAGER, MESH> Parent;
        typedef typename Parent::Base Base;
        typedef typename Parent::Wrapper Wrapper;

        typedef typename solverManagerBase::Region<Mesh> Region;
        typedef typename Parent::Control Control;
        typedef typename Parent::Settings Settings;


private:

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
            virtual void create() const = 0;


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
        virtual const Manager& manager() const = 0;

        //- Region
        virtual Region& region() const = 0;

        //- Control
        virtual Control& control() const = 0;

        //- Settings
        virtual Settings& settings() const = 0;
};


/*---------------------------------------------------------------------------*\
             Class solverManagerStorage::GlobalStorage Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER, class MESH>
class GlobalStorage
:
    public Storage<MANAGER, MESH>
{
public:

    // Public typedefs

        typedef MANAGER Manager;

        typedef typename Manager::ManagerMesh ManagerMesh;
        typedef MESH Mesh;

        typedef Storage<MANAGER, MESH> Parent;
        typedef typename Parent::Base Base;
        typedef typename Parent::Wrapper Wrapper;

        typedef typename Parent::Control Control;
        typedef typename Parent::Settings Settings;


private:

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
            virtual void create() const = 0;


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
        virtual const Manager& manager() const = 0;

        //- Control
        virtual Control& control() const = 0;

        //- Settings
        virtual Settings& settings() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solverManagerStorage

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solverManagerWrappers
{
    template <class MANAGER>
    struct RegionStorage
    {
        typedef MANAGER Manager;

        typedef typename MANAGER::ManagerMesh ManagerMesh;
        typedef typename MANAGER::Mesh Mesh;

        typedef solverManagerStorage::RegionStorage<Manager, Mesh> Type;
    };


    template <class MANAGER>
    struct GlobalStorage
    {
        typedef MANAGER Manager;

        typedef typename MANAGER::ManagerMesh ManagerMesh;
        typedef typename MANAGER::ManagerMesh Mesh;

        typedef solverManagerStorage::GlobalStorage<Manager, Mesh> Type;
    };
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
