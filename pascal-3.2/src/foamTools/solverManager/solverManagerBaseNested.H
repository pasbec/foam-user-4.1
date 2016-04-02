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
    solverManagerBaseNested

Description
    ...

\*---------------------------------------------------------------------------*/

#ifndef solverManagerBaseNested_H
#define solverManagerBaseNested_H

#include "solverManager.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define smBaseNested_new(BASEPARENT, PARENT, MANAGER, NAME)                  \
                                                                              \
    public: class NAME                                                        \
    :                                                                         \
        public solverManager<MANAGER::ManagerMesh>::BASEPARENT,               \
        public PARENT<MANAGER>                                                \
    {                                                                         \
        friend class MANAGER;                                                 \
                                                                              \
        public: typedef solverManager<MANAGER::ManagerMesh>::BASEPARENT Base; \
        public: typedef PARENT<MANAGER> BaseNested;                           \
                                                                              \
        private: NAME(const NAME&);                                           \
                                                                              \
        private: void operator=(const NAME&);                                 \
                                                                              \
        public: NAME(const MANAGER& manager)                                  \
        :                                                                     \
            Base(),                                                           \
            BaseNested(manager)                                               \
        {};                                                                   \
                                                                              \
        public: virtual ~NAME() {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define smBaseNested_end()                                                   \
                                                                              \
    };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define smBaseNested_setItem_P1(TYPE, NAME, P1)                              \
                                                                              \
                NAME##Ptr_.set(new TYPE(P1));

# define smBaseNested_setItem_P2(TYPE, NAME, P1,P2)                           \
                                                                              \
                NAME##Ptr_.set(new TYPE(P1,P2));

# define smBaseNested_setItem_P3(TYPE, NAME, P1,P2,P3)                        \
                                                                              \
                NAME##Ptr_.set(new TYPE(P1,P2,P3));

# define smBaseNested_setItem_P4(TYPE, NAME, P1,P2,P3,P4)                     \
                                                                              \
                NAME##Ptr_.set(new TYPE(P1,P2,P3,P4));

# define smBaseNested_setItem_P5(TYPE, NAME, P1,P2,P3,P4,P5)                  \
                                                                              \
                NAME##Ptr_.set(new TYPE(P1,P2,P3,P4,P5));


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define smBaseNested_addModItemNoInit_P1(TYPE, NAME, MOD, P1)                \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: MOD TYPE& NAME() const                                        \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                smBaseNested_setItem_P1(TYPE, NAME, P1)                       \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };

# define smBaseNested_addModItemNoInit_P2(TYPE, NAME, MOD, P1,P2)             \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: MOD TYPE& NAME() const                                        \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                smBaseNested_setItem_P2(TYPE, NAME, P1,P2)                    \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };

# define smBaseNested_addModItemNoInit_P3(TYPE, NAME, MOD, P1,P2,P3)          \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: MOD TYPE& NAME() const                                        \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                smBaseNested_setItem_P3(TYPE, NAME, P1,P2,P3)                 \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };

# define smBaseNested_addModItemNoInit_P4(TYPE, NAME, MOD, P1,P2,P3,P4)       \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: MOD TYPE& NAME() const                                        \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                smBaseNested_setItem_P4(TYPE, NAME, P1,P2,P3,P4)              \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };

# define smBaseNested_addModItemNoInit_P5(TYPE, NAME, MOD, P1,P2,P3,P4,P5)    \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: MOD TYPE& NAME() const                                        \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                smBaseNested_setItem_P5(TYPE, NAME, P1,P2,P3,P4,P5)           \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define smBaseNested_setInitItem_P1(TYPE, NAME, P1)                          \
                                                                              \
    smBaseNested_setItem_P1(TYPE, NAME, P1);                                  \
                                                                              \
                NAME##Ptr_->checkInit();

# define smBaseNested_setInitItem_P2(TYPE, NAME, P1,P2)                       \
                                                                              \
    smBaseNested_setItem_P2(TYPE, NAME, P1,P2);                               \
                                                                              \
                NAME##Ptr_->checkInit();

# define smBaseNested_setInitItem_P3(TYPE, NAME, P1,P2,P3)                    \
                                                                              \
    smBaseNested_setItem_P3(TYPE, NAME, P1,P2,P3);                            \
                                                                              \
                NAME##Ptr_->checkInit();

# define smBaseNested_setInitItem_P4(TYPE, NAME, P1,P2,P3,P4)                 \
                                                                              \
    smBaseNested_setItem_P4(TYPE, NAME, P1,P2,P3,P4);                         \
                                                                              \
                NAME##Ptr_->checkInit();

# define smBaseNested_setInitItem_P5(TYPE, NAME, P1,P2,P3,P4,P5)              \
                                                                              \
    smBaseNested_setItem_P5(TYPE, NAME, P1,P2,P3,P4,P5);                      \
                                                                              \
                NAME##Ptr_->checkInit();


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define smBaseNested_addModItem_P1(TYPE, NAME, MOD, P1)                      \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: MOD TYPE& NAME() const                                        \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                smBaseNested_setInitItem_P1(TYPE, NAME, P1)                   \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };

# define smBaseNested_addModItem_P2(TYPE, NAME, MOD, P1,P2)                   \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: MOD TYPE& NAME() const                                        \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                smBaseNested_setInitItem_P2(TYPE, NAME, P1,P2)                \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };

# define smBaseNested_addModItem_P3(TYPE, NAME, MOD, P1,P2,P3)                \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: MOD TYPE& NAME() const                                        \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                smBaseNested_setInitItem_P3(TYPE, NAME, P1,P2,P3)             \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };

# define smBaseNested_addModItem_P4(TYPE, NAME, MOD, P1,P2,P3,P4)             \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: MOD TYPE& NAME() const                                        \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                smBaseNested_setInitItem_P4(TYPE, NAME, P1,P2,P3,P4)          \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };

# define smBaseNested_addModItem_P5(TYPE, NAME, MOD, P1,P2,P3,P4,P5)          \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: MOD TYPE& NAME() const                                        \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                smBaseNested_setInitItem_P5(TYPE, NAME, P1,P2,P3,P4,P5)       \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define smBaseNested_addItem_P1(TYPE, NAME, P1)                              \
                                                                              \
    smBaseNested_addModItem_P1(TYPE, NAME, , P1)

# define smBaseNested_addItem_P2(TYPE, NAME, P1,P2)                           \
                                                                              \
    smBaseNested_addModItem_P2(TYPE, NAME, , P1,P2)

# define smBaseNested_addItem_P3(TYPE, NAME, P1,P2,P3)                        \
                                                                              \
    smBaseNested_addModItem_P3(TYPE, NAME, , P1,P2,P3)

# define smBaseNested_addItem_P4(TYPE, NAME, P1,P2,P3,P4)                     \
                                                                              \
    smBaseNested_addModItem_P4(TYPE, NAME, , P1,P2,P3,P4)

# define smBaseNested_addItem_P5(TYPE, NAME, P1,P2,P3,P4,P5)                  \
                                                                              \
    smBaseNested_addModItem_P5(TYPE, NAME, , P1,P2,P3,P4,P5)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define smBaseNested_addVirtualItem_P1(TYPE, NAME, P1)                       \
                                                                              \
    smBaseNested_addModItem_P1(TYPE, NAME, virtual, P1)

# define smBaseNested_addVirtualItem_P2(TYPE, NAME, P1,P2)                    \
                                                                              \
    smBaseNested_addModItem_P2(TYPE, NAME, virtual, P1,P2)

# define smBaseNested_addVirtualItem_P3(TYPE, NAME, P1,P2,P3)                 \
                                                                              \
    smBaseNested_addModItem_P3(TYPE, NAME, virtual, P1,P2,P3)

# define smBaseNested_addVirtualItem_P4(TYPE, NAME, P1,P2,P3,P4)              \
                                                                              \
    smBaseNested_addModItem_P4(TYPE, NAME, virtual, P1,P2,P3,P4)

# define smBaseNested_addVirtualItem_P5(TYPE, NAME, P1,P2,P3,P4,P5)           \
                                                                              \
    smBaseNested_addModItem_P5(TYPE, NAME, virtual, P1,P2,P3,P4,P5)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define smBaseNested_addConstItem_P1(TYPE, NAME, P1)                         \
                                                                              \
    smBaseNested_addModItem_P1(TYPE, NAME, const, P1)

# define smBaseNested_addConstItem_P2(TYPE, NAME, P1,P2)                      \
                                                                              \
    smBaseNested_addModItem_P2(TYPE, NAME, const, P1,P2)

# define smBaseNested_addConstItem_P3(TYPE, NAME, P1,P2,P3)                   \
                                                                              \
    smBaseNested_addModItem_P3(TYPE, NAME, const, P1,P2,P3)

# define smBaseNested_addConstItem_P4(TYPE, NAME, P1,P2,P3,P4)                \
                                                                              \
    smBaseNested_addModItem_P4(TYPE, NAME, const, P1,P2,P3,P4)

# define smBaseNested_addConstItem_P5(TYPE, NAME, P1,P2,P3,P4,P5)             \
                                                                              \
    smBaseNested_addModItem_P5(TYPE, NAME, const, P1,P2,P3,P4,P5)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define smBaseNested_addVirtualConstItem_P1(TYPE, NAME, P1)                  \
                                                                              \
    smBaseNested_addModItem_P1(TYPE, NAME, virtual const, P1)

# define smBaseNested_addVirtualConstItem_P2(TYPE, NAME, P1,P2)               \
                                                                              \
    smBaseNested_addModItem_P2(TYPE, NAME, virtual const, P1,P2)

# define smBaseNested_addVirtualConstItem_P3(TYPE, NAME, P1,P2,P3)            \
                                                                              \
    smBaseNested_addModItem_P3(TYPE, NAME, virtual const, P1,P2,P3)

# define smBaseNested_addVirtualConstItem_P4(TYPE, NAME, P1,P2,P3,P4)         \
                                                                              \
    smBaseNested_addModItem_P4(TYPE, NAME, virtual const, P1,P2,P3,P4)

# define smBaseNested_addVirtualConstItem_P5(TYPE, NAME, P1,P2,P3,P4,P5)      \
                                                                              \
    smBaseNested_addModItem_P5(TYPE, NAME, virtual const, P1,P2,P3,P4,P5)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define smManager_addSettings_this()                                         \
                                                                              \
    smBaseNested_addVirtualItem_P1(Settings, settings, *this)

# define smManager_addRegions_this()                                          \
                                                                              \
    smBaseNested_addVirtualItem_P1(Regions, regions, *this)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class solverManagerBaseNested Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER>
class solverManagerBaseNested
{
public:

    // Public typedefs

        typedef MANAGER Manager;
        typedef typename MANAGER::ManagerMesh ManagerMesh;
        typedef typename ManagerMesh::Mesh ThisMesh;


private:

    // Private Data

        //- Reference to args
        const argList& args_;

        //- Reference to time
        Time& time_;

        //- Reference to mesh
        ManagerMesh& mesh_;

        //- Reference to manager
        const Manager& manager_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        solverManagerBaseNested(const solverManagerBaseNested&);

        //- Disallow default bitwise assignment
        void operator=(const solverManagerBaseNested&);


protected:

    // Protected Member Functions

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
            ManagerMesh& mesh() const
            {
                return mesh_;
            }

            //- Is this the master manager?
            inline const Manager& manager() const
            {
                return manager_;
            };

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
                return manager_.propertiesDict().subDict("settings");
            };

            //- Reference to storage dictionary
            inline const dictionary& storageDict() const
            {
                return manager_.propertiesDict().subDict("storage");
            };

// TODO: Move all these basic functions in the marked block in their
// |     corresponding nested classes
// |     starting with solverManagerBaseNested. Rename
// v     solverManagerBaseNested to ?????
            //- Reference to settings dictionary for named region
            const dictionary& settingsRegionDict
            (
                const word& region = polyMesh::defaultRegion
            ) const
            {
                return manager_.settingsRegionDict(region);
            };

            //- Reference to storage dictionary for named region
            const dictionary& storageRegionDict
            (
                const word& region = polyMesh::defaultRegion
            ) const
            {
                return manager_.storageRegionDict(region);
            };

            //- Reference to dictionary of named item of storage
            bool storageItemDict
            (
                const word& item,
                dictionary& dict,
                const word& region = polyMesh::defaultRegion
            ) const
            {
                return manager_.storageItemDict(item, dict, region);
            };
// ^
// |
// |


public:

    // Constructor

        solverManagerBaseNested(const MANAGER& manager)
        :
            args_(manager.args()),
            time_(manager.time()),
            mesh_(manager.mesh()),
            manager_(manager)
        {};

    // Destructor

        virtual ~solverManagerBaseNested() {};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
