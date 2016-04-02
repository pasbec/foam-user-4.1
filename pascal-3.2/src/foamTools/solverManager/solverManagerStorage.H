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

#include "solverManager.H"
#include "solverManagerBaseNested.H"

// TODO: Macro new for regionMeshes where mesh is the mesh of one region:
//       Storage should be based on MeshType
//       with a template specification for dydnamicFvMesh

// TODO: Use mesh().name() for the virtual functions settingsRegionDict(...), storageRegionDict(...)
//       and storageItemDict(...) instead of polyMesh::defaultRegion as default
//       parameter for the region word

// TODO TODO TODO TODO :
//       Move all 'intelligent' code in class "Item"
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

# define smStorage_new(MANAGER, NAME)                                         \
                                                                              \
    smBaseNested_new(Storage, solverManagerStorage, MANAGER, NAME)            \
                                                                              \
    private: virtual void init() const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define smStorage_end()                                                      \
                                                                              \
    smBaseNested_end()


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define smStorage_add(TYPE, NAME)                                            \
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
        bool isValid = storageItemDict(name, dict);                           \
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class solverManagerStorage Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER>
class solverManagerStorage
:
    public solverManagerBaseNested<MANAGER>
{
public:

    // Public typedefs

        typedef MANAGER Manager;
        typedef typename MANAGER::ManagerMesh ManagerMesh;
        typedef typename ManagerMesh::Mesh ThisMesh;

// TODO: Items


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        solverManagerStorage(const solverManagerStorage&);

        //- Disallow default bitwise assignment
        void operator=(const solverManagerStorage&);


public:

    // Constructor

        solverManagerStorage(const MANAGER& manager)
        :
            solverManagerBaseNested<MANAGER>(manager)
        {};

    // Destructor

        virtual ~solverManagerStorage() {};
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
