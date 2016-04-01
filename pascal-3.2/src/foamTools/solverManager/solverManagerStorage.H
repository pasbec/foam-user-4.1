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

#include "autoPtr.H"
#include "argList.H"
#include "foamTime.H"
#include "IOdictionary.H"
#include "solutionControl.H"

#include "solverManager.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define solverManagerStorage_new(MANAGER, NAME)                              \
                                                                              \
    solverManagerNestedClass_new(MANAGER, Storage, NAME)                      \
                                                                              \
    protected: virtual void init(const word init = "default") const;


# define solverManagerStorage_add(TYPE, DATA)                                 \
                                                                              \
    private: mutable autoPtr<TYPE> DATA##Ptr_;                                \
                                                                              \
    private: mutable Switch DATA##Locked_;                                    \
                                                                              \
    private: mutable Switch DATA##Disabled_;                                  \
                                                                              \
    private: void make_##DATA(const word init = "default") const              \
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
                    create_##DATA(init);                                      \
                                                                              \
                    DATA##Locked_ = false;                                    \
                }                                                             \
                else                                                          \
                {                                                             \
                    FatalErrorIn                                              \
                    (manager_.type() + "Storage::make_" + #DATA + "()")       \
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
                (manager_.type() + "Storage::make_" + #DATA + "()")           \
                    << "Pointer for "                                         \
                        << #DATA                                              \
                        << " is still unallocated after creation."            \
                        << " Make sure the data member " << #DATA << "Ptr_ "  \
                        << " is set correctly in its create-function: "       \
                        << "create_" << #DATA                                 \
                        << "(const word init = 'default')!"                   \
                        << abort(FatalError);                                 \
            }                                                                 \
        }                                                                     \
        else                                                                  \
        {                                                                     \
            FatalErrorIn                                                      \
            (manager_.type() + "Storage::make_" + #DATA + "()")               \
                << "The storage item "                                        \
                    << #DATA                                                  \
                    << " is globally disabled and cannot be used."            \
                    << " Check your storage settings!"                        \
                    << abort(FatalError);                                     \
                                                                              \
        }                                                                     \
    }                                                                         \
                                                                              \
    private: void create_##DATA(const word init = "default") const;           \
                                                                              \
    public: void DATA##Enable() const                                         \
    {                                                                         \
        if (DATA##Disabled_)                                                  \
        {                                                                     \
            DATA##Disabled_ = false;                                          \
        }                                                                     \
        else                                                                  \
        {                                                                     \
            FatalErrorIn                                                      \
            (manager_.type() + "Storage::" + #DATA + "Enable()")              \
                << "The storage item "                                        \
                    << #DATA                                                  \
                    << " is already disabled!"                                \
                    << abort(FatalError);                                     \
        }                                                                     \
    };                                                                        \
                                                                              \
    public: void DATA##Disable() const                                        \
    {                                                                         \
        if (!DATA##Disabled_)                                                 \
        {                                                                     \
            if (DATA##Ptr_.valid())                                           \
            {                                                                 \
                FatalErrorIn                                                  \
                (manager_.type() + "Storage::" + #DATA + "Disable()")         \
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
    public: inline bool DATA##IsEmpty() const                                 \
    {                                                                         \
        return DATA##Ptr_.empty();                                            \
    };                                                                        \
                                                                              \
    public: inline bool DATA##IsValid() const                                 \
    {                                                                         \
        return DATA##Ptr_.valid();                                            \
    };                                                                        \
                                                                              \
    public: TYPE& DATA(const word init = "default") const                     \
    {                                                                         \
        make_##DATA(init);                                                    \
                                                                              \
        return DATA##Ptr_();                                                  \
    };                                                                        \
                                                                              \
    public: TYPE* DATA##Ptr(const word init = "default") const                \
    {                                                                         \
        make_##DATA(init);                                                    \
                                                                              \
        return DATA##Ptr_.operator->();                                       \
    };


# define solverManagerStorage_end()                                           \
                                                                              \
    solverManagerNestedClass_end()


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class solverManagerStorage Declaration
\*---------------------------------------------------------------------------*/

template <class MESH>
class solverManagerStorage
{
public:

    // Public typedefs

        typedef MESH Mesh;


private:

    // Private Data

        //- Reference to args
        const argList& args_;

        //- Reference to time
        Time& time_;

        //- Reference to mesh
        MESH& mesh_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        solverManagerStorage(const solverManagerStorage&);

        //- Disallow default bitwise assignment
        void operator=(const solverManagerStorage&);


protected:

    // Protected Member Functions

        //- Reference to args
        inline const argList& args() const
        {
            return args_;
        };

        //- Reference to time
        inline Time& time() const
        {
            return time_;
        };

        //- Return reference to mesh
        inline MESH& mesh() const
        {
            return mesh_;
        }

        //- Initialize data
        virtual void init(const word init) const = 0;


public:

    // Constructor

        //- Construct from mesh
        solverManagerStorage
        (
            const argList& args,
            Time& time,
            MESH& mesh
        )
        :
            args_(args), time_(time), mesh_(mesh)
        {};

    // Destructor

        virtual ~solverManagerStorage() {};
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
