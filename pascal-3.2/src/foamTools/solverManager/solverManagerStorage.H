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

# define solverManagerStorage_new(MANAGER, NAME)                              \
                                                                              \
    solverManagerBaseNested_new(Storage, solverManagerStorage, MANAGER, NAME) \
                                                                              \
    protected: virtual void init(const word& init = "default") const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define solverManagerStorage_end()                                           \
                                                                              \
    solverManagerBaseNested_end()


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define solverManagerStorage_add(TYPE, DATA)                                 \
                                                                              \
    private: mutable autoPtr<TYPE> DATA##Ptr_;                                \
                                                                              \
    private: mutable Switch DATA##Optional_;                                  \
                                                                              \
    private: mutable Switch DATA##Locked_;                                    \
                                                                              \
    private: mutable Switch DATA##Disabled_;                                  \
                                                                              \
    private: inline void mandatory_##DATA() const                             \
    {                                                                         \
        DATA##Optional_ = false;                                              \
    };                                                                        \
                                                                              \
    private: inline void optional_##DATA() const                              \
    {                                                                         \
        DATA##Optional_ = true;                                               \
    };                                                                        \
                                                                              \
    private: void enable_##DATA() const                                       \
    {                                                                         \
        if (DATA##Disabled_)                                                  \
        {                                                                     \
            DATA##Disabled_ = false;                                          \
        }                                                                     \
        else                                                                  \
        {                                                                     \
            FatalErrorIn                                                      \
            (this->manager().type() + "Storage::enable_"+#DATA+"()")          \
                << "The storage item "                                        \
                    << #DATA                                                  \
                    << " is already disabled!"                                \
                    << abort(FatalError);                                     \
        }                                                                     \
    };                                                                        \
                                                                              \
    private: void disable_##DATA() const                                      \
    {                                                                         \
        if (!DATA##Disabled_)                                                 \
        {                                                                     \
            if (DATA##Ptr_.valid())                                           \
            {                                                                 \
                FatalErrorIn                                                  \
                (this->manager().type() + "Storage::disable_"+#DATA+"()")     \
                    << "The storage item "                                    \
                        << #DATA                                              \
                        << " is still allocated and cannot be disabled."      \
                        << " Clear it first!"                                 \
                        << abort(FatalError);                                 \
                                                                              \
            }                                                                 \
            else                                                              \
            {                                                                 \
                DATA##Disabled_ = true;                                       \
            }                                                                 \
        }                                                                     \
    };                                                                        \
                                                                              \
    private: void create_##DATA(const word& init = "default") const           \
    {                                                                         \
        if (!DATA##Disabled_)                                                 \
        {                                                                     \
            if (DATA##Ptr_.empty())                                           \
            {                                                                 \
                if (!DATA##Locked_)                                           \
                {                                                             \
                    DATA##Locked_ = true;                                     \
                                                                              \
                    Info << "Create " << #DATA                                \
                        << " for region " << mesh().name()                    \
                        << endl;                                              \
                                                                              \
                    init_##DATA(init);                                        \
                                                                              \
                    DATA##Locked_ = false;                                    \
                }                                                             \
                else                                                          \
                {                                                             \
                    FatalErrorIn                                              \
                    (this->manager().type() + "Storage::create_"+#DATA+"()")  \
                        << "Pointer for "                                     \
                            << #DATA                                          \
                            << " is part of a circular reference."            \
                            << " Check all occurances in storage!"            \
                            << abort(FatalError);                             \
                }                                                             \
            }                                                                 \
                                                                              \
            if (DATA##Ptr_.empty())                                           \
            {                                                                 \
                FatalErrorIn                                                  \
                (this->manager().type() + "Storage::create_"+#DATA+"()")      \
                    << "Pointer for "                                         \
                        << #DATA                                              \
                        << " is still unallocated after creation."            \
                        << " Make sure the data member " << #DATA << "Ptr_ "  \
                        << " is set correctly in its init-function: "         \
                        << "init_" << #DATA << "(const word& init)!"          \
                        << abort(FatalError);                                 \
            }                                                                 \
        }                                                                     \
        else                                                                  \
        {                                                                     \
            FatalErrorIn                                                      \
            (this->manager().type() + "Storage::create_"+#DATA+"()")          \
                << "The storage item "                                        \
                    << #DATA                                                  \
                    << " is globally disabled and cannot be used."            \
                    << " Check your storage settings!"                        \
                    << abort(FatalError);                                     \
                                                                              \
        }                                                                     \
    };                                                                        \
                                                                              \
    private: bool dict_##DATA                                                 \
    (                                                                         \
        dictionary& dict,                                                     \
        const word& name = #DATA                                              \
    ) const                                                                   \
    {                                                                         \
        bool isValid = storageItemDict(name, dict);                           \
                                                                              \
        if (!DATA##Optional_)                                                 \
        {                                                                     \
            return !isValid;                                                  \
        }                                                                     \
        else                                                                  \
        {                                                                     \
            return isValid;                                                   \
        }                                                                     \
    };                                                                        \
                                                                              \
    private: bool switch_##DATA                                               \
    (                                                                         \
        dictionary& dict,                                                     \
        const word& name = #DATA                                              \
    ) const                                                                   \
    {                                                                         \
        disable_##DATA();                                                     \
                                                                              \
        if (dict_##DATA(dict, name))                                          \
        {                                                                     \
            enable_##DATA();                                                  \
                                                                              \
            return true;                                                      \
        }                                                                     \
                                                                              \
        if (!DATA##Optional_)                                                 \
        {                                                                     \
            FatalErrorIn                                                      \
            (this->manager().type() + "Storage::switch_"+#DATA+"()")          \
                        << "Storage element "                                 \
                            << #DATA                                          \
                            << " is mandatory and cannot be disabled."        \
                            << " Check settings in storage and properties!"   \
                            << abort(FatalError);                             \
        }                                                                     \
        else                                                                  \
        {                                                                     \
            Info << "Disable " << #DATA                                       \
                << " for region " << mesh().name()                            \
                << endl;                                                      \
        }                                                                     \
                                                                              \
        return false;                                                         \
    };                                                                        \
                                                                              \
    private: bool make_##DATA(const word& name = #DATA) const                 \
    {                                                                         \
        word init = "";                                                       \
                                                                              \
        dictionary dict;                                                      \
        if (switch_##DATA(dict, name))                                        \
        {                                                                     \
            init = dict.lookupOrDefault("init", word("default"));             \
                                                                              \
            create_##DATA(init);                                              \
                                                                              \
            return true;                                                      \
        }                                                                     \
                                                                              \
        return false;                                                         \
    };                                                                        \
                                                                              \
    private: bool makeMandatory_##DATA(const word& name = #DATA) const        \
    {                                                                         \
         mandatory_##DATA();                                                  \
                                                                              \
         return make_##DATA(name);                                            \
    };                                                                        \
                                                                              \
    private: bool makeOptional_##DATA(const word& name = #DATA) const         \
    {                                                                         \
         optional_##DATA();                                                   \
                                                                              \
         return make_##DATA(name);                                            \
    };                                                                        \
                                                                              \
    private: void init_##DATA(const word& init = "default") const;            \
                                                                              \
    public: inline bool is_##DATA() const                                     \
    {                                                                         \
        return !DATA##Disabled_;                                              \
    };                                                                        \
                                                                              \
    public: inline bool isNot_##DATA() const                                  \
    {                                                                         \
        return DATA##Disabled_;                                               \
    };                                                                        \
                                                                              \
    public: inline bool isEmpty_##DATA() const                                \
    {                                                                         \
        return DATA##Ptr_.empty();                                            \
    };                                                                        \
                                                                              \
    public: inline bool isValid_##DATA() const                                \
    {                                                                         \
        return DATA##Ptr_.valid();                                            \
    };                                                                        \
                                                                              \
    public: TYPE& DATA(const word& init = "default") const                    \
    {                                                                         \
        create_##DATA(init);                                                  \
                                                                              \
        return DATA##Ptr_();                                                  \
    };                                                                        \
                                                                              \
    public: TYPE* DATA##Ptr(const word& init = "default") const               \
    {                                                                         \
        create_##DATA(init);                                                  \
                                                                              \
        return DATA##Ptr_.operator->();                                       \
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
