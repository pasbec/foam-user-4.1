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

// TODO: makeOptional_ has to be connected to settings!

// TODO: Implement Storage Items
//
//       Move all 'more or less intelligent' code in class "Item"
//       and use only macro to get all items in
//       the Storage namespace e.g.:
//
//         autoPtr<Item> pItemPtr_;
//
//       And then take the macros only for:
//
//         p() => pItemPtr->get(),
//         pPtr() => pItemPtr_->Ptr(),
//         is_p() => pItemPtr.is(), ...

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGIONSTORAGE_NEW(MANAGER)                                        \
                                                                              \
    SM_WRAPPER_NEW(RegionStorage, MANAGER, Storage)                           \
                                                                              \
    private: virtual void init() const;


# define SM_GLOBALSTORAGE_NEW(MANAGER)                                        \
                                                                              \
    SM_WRAPPER_NEW(GlobalStorage, MANAGER, Storage)                           \
                                                                              \
    private: virtual void init() const;


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

// TODO: Use Base class for itemDict?

/*---------------------------------------------------------------------------*\
              Class solverManagerStorage::RegionStorage Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER, class MESH>
class RegionStorage
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


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        RegionStorage(const RegionStorage&);

        //- Disallow default bitwise assignment
        void operator=(const RegionStorage&);

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

        RegionStorage(const Manager& manager, Mesh& mesh)
        :
            Base(),
            Wrapper(manager, mesh)
        {};

    // Destructor

        virtual ~RegionStorage() {};

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
                return storageDict();
            };


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
             Class solverManagerStorage::GlobalStorage Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER, class MESH>
class GlobalStorage
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


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        GlobalStorage(const GlobalStorage&);

        //- Disallow default bitwise assignment
        void operator=(const GlobalStorage&);

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

        GlobalStorage(const Manager& manager, Mesh& mesh)
        :
            Base(),
            Wrapper(manager, mesh)
        {};

    // Destructor

        virtual ~GlobalStorage() {};

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
                return storageDict();
            };


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
