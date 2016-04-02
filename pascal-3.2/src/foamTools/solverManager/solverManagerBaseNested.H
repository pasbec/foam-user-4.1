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

# define solverManagerBaseNested_new(BASEPARENT, PARENT, MANAGER, NAME)       \
                                                                              \
    public: class NAME                                                        \
    :                                                                         \
        public solverManager<MANAGER::ManagerMesh>::BASEPARENT,               \
        public PARENT<MANAGER>                                                \
    {                                                                         \
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
        {                                                                     \
            init();                                                           \
        };                                                                    \
                                                                              \
        public: virtual ~NAME() {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define solverManagerBaseNested_end()                                        \
                                                                              \
    };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define solverManagerBaseNested_setItemP1(TYPE, ITEM, P1)                    \
                                                                              \
                ITEM##Ptr_.set(new TYPE(P1));

# define solverManagerBaseNested_setItemP2(TYPE, ITEM, P1,P2)                 \
                                                                              \
                ITEM##Ptr_.set(new TYPE(P1,P2));

# define solverManagerBaseNested_setItemP3(TYPE, ITEM, P1,P2,P3)              \
                                                                              \
                ITEM##Ptr_.set(new TYPE(P1,P2,P3));

# define solverManagerBaseNested_setItemP4(TYPE, ITEM, P1,P2,P3,P4)           \
                                                                              \
                ITEM##Ptr_.set(new TYPE(P1,P2,P3,P4));

# define solverManagerBaseNested_setItemP5(TYPE, ITEM, P1,P2,P3,P4,P5)        \
                                                                              \
                ITEM##Ptr_.set(new TYPE(P1,P2,P3,P4,P5));


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define solverManagerBaseNested_addItemModP1(TYPE, ITEM, MOD, P1)            \
                                                                              \
        private: mutable autoPtr<TYPE> ITEM##Ptr_;                            \
                                                                              \
        public: MOD TYPE& ITEM() const                                        \
        {                                                                     \
            if (ITEM##Ptr_.empty())                                           \
            {                                                                 \
                solverManagerBaseNested_setItemP1(TYPE, ITEM, P1);            \
            }                                                                 \
                                                                              \
            return ITEM##Ptr_();                                              \
        }

# define solverManagerBaseNested_addItemModP2(TYPE, ITEM, MOD, P1,P2)         \
                                                                              \
        private: mutable autoPtr<TYPE> ITEM##Ptr_;                            \
                                                                              \
        public: MOD TYPE& ITEM() const                                        \
        {                                                                     \
            if (ITEM##Ptr_.empty())                                           \
            {                                                                 \
                solverManagerBaseNested_setItemP2(TYPE, ITEM, P1,P2);         \
            }                                                                 \
                                                                              \
            return ITEM##Ptr_();                                              \
        }

# define solverManagerBaseNested_addItemModP3(TYPE, ITEM, MOD, P1,P2,P3)      \
                                                                              \
        private: mutable autoPtr<TYPE> ITEM##Ptr_;                            \
                                                                              \
        public: MOD TYPE& ITEM() const                                        \
        {                                                                     \
            if (ITEM##Ptr_.empty())                                           \
            {                                                                 \
                solverManagerBaseNested_setItemP3(TYPE, ITEM, P1,P2,P3);      \
            }                                                                 \
                                                                              \
            return ITEM##Ptr_();                                              \
        }

# define solverManagerBaseNested_addItemModP4(TYPE, ITEM, MOD, P1,P2,P3,P4)   \
                                                                              \
        private: mutable autoPtr<TYPE> ITEM##Ptr_;                            \
                                                                              \
        public: MOD TYPE& ITEM() const                                        \
        {                                                                     \
            if (ITEM##Ptr_.empty())                                           \
            {                                                                 \
                solverManagerBaseNested_setItemP4(TYPE, ITEM, P1,P2,P3,P4);   \
            }                                                                 \
                                                                              \
            return ITEM##Ptr_();                                              \
        }

# define solverManagerBaseNested_addItemModP5(TYPE, ITEM, MOD, P1,P2,P3,P4,P5)\
                                                                              \
        private: mutable autoPtr<TYPE> ITEM##Ptr_;                            \
                                                                              \
        public: MOD TYPE& ITEM() const                                        \
        {                                                                     \
            if (ITEM##Ptr_.empty())                                           \
            {                                                                 \
                solverManagerBaseNested_setItemP5(TYPE, ITEM, P1,P2,P3,P4,P5);\
            }                                                                 \
                                                                              \
            return ITEM##Ptr_();                                              \
        }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define solverManagerBaseNested_addItemP1(TYPE, ITEM, P1)                    \
                                                                              \
        solverManagerBaseNested_addItemModP1                                  \
        (TYPE, ITEM, , P1);

# define solverManagerBaseNested_addItemP2(TYPE, ITEM, P1,P2)                 \
                                                                              \
        solverManagerBaseNested_addItemModP2                                  \
        (TYPE, ITEM, , P1,P2);

# define solverManagerBaseNested_addItemP3(TYPE, ITEM, P1,P2,P3)              \
                                                                              \
        solverManagerBaseNested_addItemModP3                                  \
        (TYPE, ITEM, , P1,P2,P3);

# define solverManagerBaseNested_addItemP4(TYPE, ITEM, P1,P2,P3,P4)           \
                                                                              \
        solverManagerBaseNested_addItemModP4                                  \
        (TYPE, ITEM, , P1,P2,P3,P4);

# define solverManagerBaseNested_addItemP5(TYPE, ITEM, P1,P2,P3,P4,P5)        \
                                                                              \
        solverManagerBaseNested_addItemModP5                                  \
        (TYPE, ITEM, , P1,P2,P3,P4,P5);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define solverManagerBaseNested_addVirtualItemP1(TYPE, ITEM, P1)             \
                                                                              \
        solverManagerBaseNested_addItemModP1                                  \
        (TYPE, ITEM, virtual, P1);

# define solverManagerBaseNested_addVirtualItemP2(TYPE, ITEM, P1,P2)          \
                                                                              \
        solverManagerBaseNested_addItemModP2                                  \
        (TYPE, ITEM, virtual, P1,P2);

# define solverManagerBaseNested_addVirtualItemP3(TYPE, ITEM, P1,P2,P3)       \
                                                                              \
        solverManagerBaseNested_addItemModP3                                  \
        (TYPE, ITEM, virtual, P1,P2,P3);

# define solverManagerBaseNested_addVirtualItemP4(TYPE, ITEM, P1,P2,P3,P4)    \
                                                                              \
        solverManagerBaseNested_addItemModP4                                  \
        (TYPE, ITEM, virtual, P1,P2,P3,P4);

# define solverManagerBaseNested_addVirtualItemP5(TYPE, ITEM, P1,P2,P3,P4,P5) \
                                                                              \
        solverManagerBaseNested_addItemModP5                                  \
        (TYPE, ITEM, virtual, P1,P2,P3,P4,P5);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define solverManagerBaseNested_addConstItemP1(TYPE, ITEM, P1)               \
                                                                              \
        solverManagerBaseNested_addItemModP1                                  \
        (TYPE, ITEM, const, P1);

# define solverManagerBaseNested_addConstItemP2(TYPE, ITEM, P1,P2)            \
                                                                              \
        solverManagerBaseNested_addItemModP2                                  \
        (TYPE, ITEM, const, P1,P2);

# define solverManagerBaseNested_addConstItemP3(TYPE, ITEM, P1,P2,P3)         \
                                                                              \
        solverManagerBaseNested_addItemModP3                                  \
        (TYPE, ITEM, const, P1,P2,P3);

# define solverManagerBaseNested_addConstItemP4(TYPE, ITEM, P1,P2,P3,P4)      \
                                                                              \
        solverManagerBaseNested_addItemModP4                                  \
        (TYPE, ITEM, const, P1,P2,P3,P4);

# define solverManagerBaseNested_addConstItemP5(TYPE, ITEM, P1,P2,P3,P4,P5)   \
                                                                              \
        solverManagerBaseNested_addItemModP5                                  \
        (TYPE, ITEM, const, P1,P2,P3,P4,P5);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define solverManagerBaseNested_addVirtualConstItemP1                        \
(TYPE, ITEM, P1)                                                              \
                                                                              \
        solverManagerBaseNested_addItemModP1                                  \
        (TYPE, ITEM, virtual const, P1);

# define solverManagerBaseNested_addVirtualConstItemP2                        \
(TYPE, ITEM, P1,P2)                                                           \
                                                                              \
        solverManagerBaseNested_addItemModP2                                  \
        (TYPE, ITEM, virtual const, P1,P2);

# define solverManagerBaseNested_addVirtualConstItemP3                        \
(TYPE, ITEM, P1,P2,P3)                                                        \
                                                                              \
        solverManagerBaseNested_addItemModP3                                  \
        (TYPE, ITEM, virtual const, P1,P2,P3);

# define solverManagerBaseNested_addVirtualConstItemP4                        \
(TYPE, ITEM, P1,P2,P3,P4)                                                     \
                                                                              \
        solverManagerBaseNested_addItemModP4                                  \
        (TYPE, ITEM, virtual const, P1,P2,P3,P4);

# define solverManagerBaseNested_addVirtualConstItemP5                        \
(TYPE, ITEM, P1,P2,P3,P4,P5)                                                  \
                                                                              \
        solverManagerBaseNested_addItemModP5                                  \
        (TYPE, ITEM, virtual const, P1,P2,P3,P4,P5);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


# define solverManagerBaseNested_addPmanager(TYPE, ITEM)                      \
                                                                              \
        solverManagerBaseNested_addItemP1                                     \
        (TYPE, ITEM, this->manager())

# define solverManagerBaseNested_addVirtualPmanager(TYPE, ITEM)               \
                                                                              \
        solverManagerBaseNested_addVirtualItemP1                              \
        (TYPE, ITEM, this->manager())

# define solverManagerBaseNested_addConstPmanager(TYPE, ITEM)                 \
                                                                              \
        solverManagerBaseNested_addConstItemP1                                \
        (TYPE, ITEM, this->manager())

# define solverManagerBaseNested_addVirtualConstPmanager(TYPE, ITEM)          \
                                                                              \
        solverManagerBaseNested_addVirtualConstItemP1                         \
        (TYPE, ITEM, this->manager())


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

            //- Reference to settings dictionary for named region
            const dictionary& settingsRegionDict
            (
                const word& region = polyMesh::defaultRegion
            ) const
            {
                return manager_.settingsRegionDict(region);
            };

            //- Reference to storage dictionary
            inline const dictionary& storageDict() const
            {
                return manager_.propertiesDict().subDict("storage");
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

        // Evolution

            //- Init sequence
            virtual void init(const word& init = "default") const {};


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
