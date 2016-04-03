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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_WRAPPER_NEW(BASE, MANAGER, NAME)                                  \
                                                                              \
    public: class NAME                                                        \
    :                                                                         \
        public solverManager##BASE<MANAGER>,                                  \
        public solverManager<ManagerMesh>::BASE                               \
    {                                                                         \
        friend class MANAGER;                                                 \
                                                                              \
        public: typedef MANAGER Manager;                                      \
        public: typedef solverManager##BASE<MANAGER> Wrapper;                 \
                                                                              \
        public: typedef typename solverManager<ManagerMesh>::BASE Base;       \
        public: typedef typename Wrapper::ManagerMesh ManagerMesh;            \
        public: typedef typename Wrapper::Mesh Mesh;                          \
                                                                              \
        private: NAME(const NAME&);                                           \
        private: void operator=(const NAME&);                                 \
                                                                              \
        public: NAME                                                          \
        (                                                                     \
            const Manager& manager,                                           \
            Mesh& mesh                                                        \
        )                                                                     \
        :                                                                     \
            Wrapper(manager, mesh),                                           \
            Base()                                                            \
        {};                                                                   \
                                                                              \
        public: virtual ~NAME() {};


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

# define SM_ADD_GLOBALSETTINGS()                                              \
                                                                              \
    SM_WRAPPER_MEMBER_VIRTUAL_ADD_P2(                                         \
        GlobalSettings, globalSettings, *this, this->mesh())

# define SM_ADD_REGIONS()                                                     \
                                                                              \
    SM_WRAPPER_MEMBER_VIRTUAL_ADD_P2(                                         \
        Regions, regions, *this, this->mesh())                                \
                                                                              \
        typedef typename DefaultRegion::Control Control;                      \
                                                                              \
        public: inline DefaultRegion::Control& control() const                \
        {                                                                     \
            return regions().defaultRegion().control();                       \
        };                                                                    \
                                                                              \
        typedef typename DefaultRegion::Settings Settings;                    \
                                                                              \
        public: inline DefaultRegion::Settings& settings() const              \
        {                                                                     \
            return regions().defaultRegion().settings();                      \
        };                                                                    \
                                                                              \
        typedef typename DefaultRegion::Storage Storage;                      \
                                                                              \
        public: inline DefaultRegion::Storage& storage() const                \
        {                                                                     \
            return regions().defaultRegion().storage();                       \
        };


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

        typedef typename MANAGER::ManagerMesh ManagerMesh;
        typedef MESH Mesh;

        typedef solverManagerWrapper Wrapper;


private:

    // Private Data

        //- Reference to manager
        const Manager& manager_;

        // Basic

            //- Reference to args
            const argList& args_;

            //- Reference to time
            Time& time_;

            //- Reference to mesh
            Mesh& mesh_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        solverManagerWrapper(const solverManagerWrapper&);

        //- Disallow default bitwise assignment
        void operator=(const solverManagerWrapper&);


public:

    // Constructor

        solverManagerWrapper
        (
            const Manager& manager,
            Mesh& mesh
        )
        :
            manager_(manager),
            args_(manager.args()),
            time_(manager.time()),
            mesh_(manager.mesh())
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
            inline Mesh& mesh() const
            {
                return mesh_;
            }

            //- Is this the master manager?
            inline const bool& master() const
            {
                return manager_.master();
            };

            //- Refrence to control dictionary
            inline const dictionary& controlDict() const
            {
                return time_.controlDict();
            }

            //- Reference to properties dictionary
            inline const dictionary& propertiesDict() const
            {
                return manager_.propertiesDict();
            };

            //- Reference to settings dictionary
            inline const dictionary& settingsDict() const
            {
                return propertiesDict().subDict("settings");
            };

            //- Reference to storage dictionary
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
