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
    solverManagerWrapper

Description
    ...

\*---------------------------------------------------------------------------*/

#ifndef solverManagerWrapper_H
#define solverManagerWrapper_H

#include "argList.H"
#include "foamTime.H"
#include "dictionary.H"

#include "solverManagerCore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_WRAPPER_NEWHEADER(WRAPPER, BASE, MANAGER, NAME)                   \
                                                                              \
    public: class NAME                                                        \
    :                                                                         \
        public solverManager##WRAPPER<MANAGER>                                \
    {                                                                         \
        friend class MANAGER;                                                 \
                                                                              \
        public: typedef MANAGER Manager;                                      \
        public: typedef solverManager##WRAPPER<MANAGER> Wrapper;              \
                                                                              \
        public: typedef typename solverManagerCore<ManagerMesh>::BASE Base;   \
        public: typedef typename Wrapper::ManagerMesh ManagerMesh;            \
        public: typedef typename Wrapper::Mesh Mesh;                          \
                                                                              \
        private: NAME(const NAME&);                                           \
        private: void operator=(const NAME&);


# define SM_WRAPPER_NEWCONSTRUCTOR(WRAPPER, BASE, MANAGER, NAME)              \
                                                                              \
        public: NAME                                                          \
        (                                                                     \
            const MANAGER& manager,                                           \
            solverManager##WRAPPER<MANAGER>::Mesh& mesh                       \
        )                                                                     \
        :                                                                     \
            solverManager##WRAPPER<MANAGER>(manager, mesh)                    \
        {};


# define SM_WRAPPER_NEWDESTRUCTOR(NAME)                                       \
                                                                              \
        public: virtual ~NAME() {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_WRAPPER_NEW(WRAPPER, BASE, MANAGER, NAME)                         \
                                                                              \
    SM_WRAPPER_NEWHEADER(WRAPPER, BASE, MANAGER, NAME)                        \
    SM_WRAPPER_NEWCONSTRUCTOR(WRAPPER, BASE, MANAGER, NAME)                   \
    SM_WRAPPER_NEWDESTRUCTOR(NAME)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_WRAPPER_END()                                                     \
                                                                              \
    };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_WRAPPER_MEMBER_SET_P1(TYPE, NAME, P1)                             \
                                                                              \
                NAME##Ptr_.set(new TYPE(P1));

# define SM_WRAPPER_MEMBER_SET_P2(TYPE, NAME, P1,P2)                          \
                                                                              \
                NAME##Ptr_.set(new TYPE(P1,P2));

# define SM_WRAPPER_MEMBER_SET_P3(TYPE, NAME, P1,P2,P3)                       \
                                                                              \
                NAME##Ptr_.set(new TYPE(P1,P2,P3));

# define SM_WRAPPER_MEMBER_SET_P4(TYPE, NAME, P1,P2,P3,P4)                    \
                                                                              \
                NAME##Ptr_.set(new TYPE(P1,P2,P3,P4));

# define SM_WRAPPER_MEMBER_SET_P5(TYPE, NAME, P1,P2,P3,P4,P5)                 \
                                                                              \
                NAME##Ptr_.set(new TYPE(P1,P2,P3,P4,P5));


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_WRAPPER_MEMBER_MOD_NOINIT_ADD_P1(TYPE, NAME, MOD, P1)             \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: MOD TYPE& NAME() const                                        \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                SM_WRAPPER_MEMBER_SET_P1(TYPE, NAME, P1)                      \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };

# define SM_WRAPPER_MEMBER_MOD_NOINIT_ADD_P2(TYPE, NAME, MOD, P1,P2)          \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: MOD TYPE& NAME() const                                        \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                SM_WRAPPER_MEMBER_SET_P2(TYPE, NAME, P1,P2)                   \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };

# define SM_WRAPPER_MEMBER_MOD_NOINIT_ADD_P3(TYPE, NAME, MOD, P1,P2,P3)       \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: MOD TYPE& NAME() const                                        \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                SM_WRAPPER_MEMBER_SET_P3(TYPE, NAME, P1,P2,P3)                \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };

# define SM_WRAPPER_MEMBER_MOD_NOINIT_ADD_P4(TYPE, NAME, MOD, P1,P2,P3,P4)    \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: MOD TYPE& NAME() const                                        \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                SM_WRAPPER_MEMBER_SET_P4(TYPE, NAME, P1,P2,P3,P4)             \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };

# define SM_WRAPPER_MEMBER_MOD_NOINIT_ADD_P5(TYPE, NAME, MOD, P1,P2,P3,P4,P5) \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: MOD TYPE& NAME() const                                        \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                SM_WRAPPER_MEMBER_SET_P5(TYPE, NAME, P1,P2,P3,P4,P5)          \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_WRAPPER_MEMBER_SETINIT_P1(TYPE, NAME, P1)                         \
                                                                              \
    SM_WRAPPER_MEMBER_SET_P1(TYPE, NAME, P1);                                 \
                                                                              \
                NAME##Ptr_->checkInit();

# define SM_WRAPPER_MEMBER_SETINIT_P2(TYPE, NAME, P1,P2)                      \
                                                                              \
    SM_WRAPPER_MEMBER_SET_P2(TYPE, NAME, P1,P2);                              \
                                                                              \
                NAME##Ptr_->checkInit();

# define SM_WRAPPER_MEMBER_SETINIT_P3(TYPE, NAME, P1,P2,P3)                   \
                                                                              \
    SM_WRAPPER_MEMBER_SET_P3(TYPE, NAME, P1,P2,P3);                           \
                                                                              \
                NAME##Ptr_->checkInit();

# define SM_WRAPPER_MEMBER_SETINIT_P4(TYPE, NAME, P1,P2,P3,P4)                \
                                                                              \
    SM_WRAPPER_MEMBER_SET_P4(TYPE, NAME, P1,P2,P3,P4);                        \
                                                                              \
                NAME##Ptr_->checkInit();

# define SM_WRAPPER_MEMBER_SETINIT_P5(TYPE, NAME, P1,P2,P3,P4,P5)             \
                                                                              \
    SM_WRAPPER_MEMBER_SET_P5(TYPE, NAME, P1,P2,P3,P4,P5);                     \
                                                                              \
                NAME##Ptr_->checkInit();


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_WRAPPER_MEMBER_MOD_ADD_P1(TYPE, NAME, MOD, P1)                    \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: MOD TYPE& NAME() const                                        \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                SM_WRAPPER_MEMBER_SETINIT_P1(TYPE, NAME, P1)                  \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };

# define SM_WRAPPER_MEMBER_MOD_ADD_P2(TYPE, NAME, MOD, P1,P2)                 \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: MOD TYPE& NAME() const                                        \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                SM_WRAPPER_MEMBER_SETINIT_P2(TYPE, NAME, P1,P2)               \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };

# define SM_WRAPPER_MEMBER_MOD_ADD_P3(TYPE, NAME, MOD, P1,P2,P3)              \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: MOD TYPE& NAME() const                                        \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                SM_WRAPPER_MEMBER_SETINIT_P3(TYPE, NAME, P1,P2,P3)            \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };

# define SM_WRAPPER_MEMBER_MOD_ADD_P4(TYPE, NAME, MOD, P1,P2,P3,P4)           \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: MOD TYPE& NAME() const                                        \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                SM_WRAPPER_MEMBER_SETINIT_P4(TYPE, NAME, P1,P2,P3,P4)         \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };

# define SM_WRAPPER_MEMBER_MOD_ADD_P5(TYPE, NAME, MOD, P1,P2,P3,P4,P5)        \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: MOD TYPE& NAME() const                                        \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                SM_WRAPPER_MEMBER_SETINIT_P5(TYPE, NAME, P1,P2,P3,P4,P5)        \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_WRAPPER_MEMBER_ADD_P1(TYPE, NAME, P1)                             \
                                                                              \
    SM_WRAPPER_MEMBER_MOD_ADD_P1(TYPE, NAME, , P1)

# define SM_WRAPPER_MEMBER_ADD_P2(TYPE, NAME, P1,P2)                          \
                                                                              \
    SM_WRAPPER_MEMBER_MOD_ADD_P2(TYPE, NAME, , P1,P2)

# define SM_WRAPPER_MEMBER_ADD_P3(TYPE, NAME, P1,P2,P3)                       \
                                                                              \
    SM_WRAPPER_MEMBER_MOD_ADD_P3(TYPE, NAME, , P1,P2,P3)

# define SM_WRAPPER_MEMBER_ADD_P4(TYPE, NAME, P1,P2,P3,P4)                    \
                                                                              \
    SM_WRAPPER_MEMBER_MOD_ADD_P4(TYPE, NAME, , P1,P2,P3,P4)

# define SM_WRAPPER_MEMBER_ADD_P5(TYPE, NAME, P1,P2,P3,P4,P5)                 \
                                                                              \
    SM_WRAPPER_MEMBER_MOD_ADD_P5(TYPE, NAME, , P1,P2,P3,P4,P5)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_WRAPPER_MEMBER_VIRTUAL_ADD_P1(TYPE, NAME, P1)                     \
                                                                              \
    SM_WRAPPER_MEMBER_MOD_ADD_P1(TYPE, NAME, virtual, P1)

# define SM_WRAPPER_MEMBER_VIRTUAL_ADD_P2(TYPE, NAME, P1,P2)                  \
                                                                              \
    SM_WRAPPER_MEMBER_MOD_ADD_P2(TYPE, NAME, virtual, P1,P2)

# define SM_WRAPPER_MEMBER_VIRTUAL_ADD_P3(TYPE, NAME, P1,P2,P3)               \
                                                                              \
    SM_WRAPPER_MEMBER_MOD_ADD_P3(TYPE, NAME, virtual, P1,P2,P3)

# define SM_WRAPPER_MEMBER_VIRTUAL_ADD_P4(TYPE, NAME, P1,P2,P3,P4)            \
                                                                              \
    SM_WRAPPER_MEMBER_MOD_ADD_P4(TYPE, NAME, virtual, P1,P2,P3,P4)

# define SM_WRAPPER_MEMBER_VIRTUAL_ADD_P5(TYPE, NAME, P1,P2,P3,P4,P5)         \
                                                                              \
    SM_WRAPPER_MEMBER_MOD_ADD_P5(TYPE, NAME, virtual, P1,P2,P3,P4,P5)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_WRAPPER_MEMBER_CONST_ADD_P1(TYPE, NAME, P1)                       \
                                                                              \
    SM_WRAPPER_MEMBER_MOD_ADD_P1(TYPE, NAME, const, P1)

# define SM_WRAPPER_MEMBER_CONST_ADD_P2(TYPE, NAME, P1,P2)                    \
                                                                              \
    SM_WRAPPER_MEMBER_MOD_ADD_P2(TYPE, NAME, const, P1,P2)

# define SM_WRAPPER_MEMBER_CONST_ADD_P3(TYPE, NAME, P1,P2,P3)                 \
                                                                              \
    SM_WRAPPER_MEMBER_MOD_ADD_P3(TYPE, NAME, const, P1,P2,P3)

# define SM_WRAPPER_MEMBER_CONST_ADD_P4(TYPE, NAME, P1,P2,P3,P4)              \
                                                                              \
    SM_WRAPPER_MEMBER_MOD_ADD_P4(TYPE, NAME, const, P1,P2,P3,P4)

# define SM_WRAPPER_MEMBER_CONST_ADD_P5(TYPE, NAME, P1,P2,P3,P4,P5)           \
                                                                              \
    SM_WRAPPER_MEMBER_MOD_ADD_P5(TYPE, NAME, const, P1,P2,P3,P4,P5)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define sSM_WRAPPER_MEMBER_VIRTUALCONST_ADD_P1(TYPE, NAME, P1)               \
                                                                              \
    SM_WRAPPER_MEMBER_MOD_ADD_P1(TYPE, NAME, virtual const, P1)

# define sSM_WRAPPER_MEMBER_VIRTUALCONST_ADD_P2(TYPE, NAME, P1,P2)            \
                                                                              \
    SM_WRAPPER_MEMBER_MOD_ADD_P2(TYPE, NAME, virtual const, P1,P2)

# define sSM_WRAPPER_MEMBER_VIRTUALCONST_ADD_P3(TYPE, NAME, P1,P2,P3)         \
                                                                              \
    SM_WRAPPER_MEMBER_MOD_ADD_P3(TYPE, NAME, virtual const, P1,P2,P3)

# define sSM_WRAPPER_MEMBER_VIRTUALCONST_ADD_P4(TYPE, NAME, P1,P2,P3,P4)      \
                                                                              \
    SM_WRAPPER_MEMBER_MOD_ADD_P4(TYPE, NAME, virtual const, P1,P2,P3,P4)

# define sSM_WRAPPER_MEMBER_VIRTUALCONST_ADD_P5(TYPE, NAME, P1,P2,P3,P4,P5)   \
                                                                              \
    SM_WRAPPER_MEMBER_MOD_ADD_P5(TYPE, NAME, virtual const, P1,P2,P3,P4,P5)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_ADD_GLOBALCONTROL_P1(TYPE, P1)                                    \
                                                                              \
        public: typedef TYPE Control;                                         \
                                                                              \
        SM_WRAPPER_MEMBER_MOD_NOINIT_ADD_P1(                                  \
            Control, control, virtual, P1)

# define SM_ADD_GLOBALCONTROL_P2(TYPE, P1,P2)                                 \
                                                                              \
        public: typedef TYPE Control;                                         \
                                                                              \
        SM_WRAPPER_MEMBER_MOD_NOINIT_ADD_P2(                                  \
            Control, control, virtual, P1,P2)

# define SM_ADD_GLOBALCONTROL_P3(TYPE, P1,P2,P3)                              \
                                                                              \
        public: typedef TYPE Control;                                         \
                                                                              \
        SM_WRAPPER_MEMBER_MOD_NOINIT_ADD_P3(                                  \
            Control, control, virtual, P1,P2,P3)

# define SM_ADD_GLOBALCONTROL_P4(TYPE, P1,P2,P3,P4)                           \
                                                                              \
        public: typedef TYPE Control;                                         \
                                                                              \
        SM_WRAPPER_MEMBER_MOD_NOINIT_ADD_P4(                                  \
            Control, control, virtual, P1,P2,P3,P4)

# define SM_ADD_GLOBALCONTROL_P5(TYPE, P1,P2,P3,P4,P5)                        \
                                                                              \
        public: typedef TYPE Control;                                         \
                                                                              \
        SM_WRAPPER_MEMBER_MOD_NOINIT_ADD_P5(                                  \
            Control, control, virtual, P1,P2,P3,P4,P5)


# define SM_GLOBALCONTROL_ADD_P1(TYPE, P1)                                    \
                                                                              \
        SM_ADD_GLOBALCONTROL_P1(TYPE, P1)

# define SM_GLOBALCONTROL_ADD_P2(TYPE, P1,P2)                                 \
                                                                              \
        SM_ADD_GLOBALCONTROL_P2(TYPE, P1,P2)

# define SM_GLOBALCONTROL_ADD_P3(TYPE, P1,P2,P3)                              \
                                                                              \
        SM_ADD_GLOBALCONTROL_P3(TYPE, P1,P2,P3)

# define SM_GLOBALCONTROL_ADD_P4(TYPE, P1,P2,P3,P4)                           \
                                                                              \
        SM_ADD_GLOBALCONTROL_P4(TYPE, P1,P2,P3,P4)

# define SM_GLOBALCONTROL_ADD_P5(TYPE, P1,P2,P3,P4,P5)                        \
                                                                              \
        SM_ADD_GLOBALCONTROL_P5(TYPE, P1,P2,P3,P4,P5)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_ADD_GLOBALSETTINGS()                                              \
                                                                              \
    SM_WRAPPER_MEMBER_VIRTUAL_ADD_P2(                                         \
        Settings, settings, *this, this->mesh())


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_ADD_GLOBALSTORAGE()                                               \
                                                                              \
        SM_WRAPPER_MEMBER_VIRTUAL_ADD_P2(                                     \
            Storage, storage, *this, this->mesh())


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_ADD_REGIONS()                                                     \
                                                                              \
    SM_WRAPPER_MEMBER_VIRTUAL_ADD_P2(                                         \
        Regions, regions, *this, this->mesh())


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class solverManagerWrapper Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER, class MESH>
class solverManagerWrapper
{
public:

    // Public typedefs

        typedef MANAGER Manager;

        typedef typename Manager::ManagerMesh ManagerMesh;
        typedef MESH Mesh;

        typedef solverManagerWrapper Wrapper;


private:

    // Private Data

        //- Reference to manager
        const Manager& manager_;

        // Basic

            //- Reference to mesh
            Mesh& mesh_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        solverManagerWrapper(const solverManagerWrapper&);

        //- Disallow default bitwise assignment
        void operator=(const solverManagerWrapper&);


protected:

    // Protected Member Functions

        // Basic

            //- Mesh
            inline Mesh& wrappedMesh() const
            {
                return mesh_;
            }


public:

    // Constructor

        solverManagerWrapper
        (
            const Manager& manager,
            Mesh& mesh
        )
        :
            manager_(manager),
            mesh_(mesh)
        {};

    // Destructor

        virtual ~solverManagerWrapper() {};

    // Member Functions

        //- Reference to manager
        inline const Manager& manager() const
        {
            return manager_;
        };

        // Basic

            //- Main arguments
            inline const argList& args() const
            {
                return manager_.args();
            };

            //- Time
            inline Time& time() const
            {
                return manager_.time();
            };

            //- Is this the master manager?
            inline const bool& master() const
            {
                return manager_.master();
            };

        // Dictionaries

            //- Control dictionary
            inline const dictionary& controlDict() const
            {
                return time().controlDict();
            }

            //- Properties dictionary
            inline const dictionary& propertiesDict() const
            {
                return manager_.propertiesDict();
            };

            //- Settings dictionary
            inline const dictionary& settingsDict() const
            {
                return propertiesDict().subDict("settings");
            };

            //- Storage dictionary
            inline const dictionary& storageDict() const
            {
                return propertiesDict().subDict("storage");
            };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
