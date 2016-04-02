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
    solverManager

Description
    ...

SourceFiles
    solverManager.C

\*---------------------------------------------------------------------------*/

#ifndef solverManager_H
#define solverManager_H

#include "solverManagerMessages.H"
#include "solverManagerRegions.H"
#include "solverManagerStorage.H"

#include "autoPtr.H"
#include "argList.H"
#include "foamTime.H"
#include "IOdictionary.H"
#include "solutionControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// TODO: Remake these macros

# define createManager(NAME, MESH)                                            \
                                                                              \
    NAME##Manager manager(MESH);


# define setManagerScope(NAME)                                                \
                                                                              \
    NAME##Manager& NAME = manager.NAME();                                     \
    NAME##Manager& manager = NAME;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define solverManagerNestedClass_new(MANAGER, TYPE, NAME)                    \
                                                                              \
    public: class NAME                                                        \
    :                                                                         \
        public solverManager<MANAGER::Mesh>::TYPE                             \
    {                                                                         \
        private: NAME(const NAME&);                                           \
                                                                              \
        private: void operator=(const NAME&);                                 \
                                                                              \
        private: const MANAGER& manager_;                                     \
                                                                              \
        protected: inline const MANAGER& manager() const                      \
        {                                                                     \
            return manager_;                                                  \
        };                                                                    \
                                                                              \
        protected: inline const dictionary& propertiesDict() const            \
        {                                                                     \
            return manager_.propertiesDict();                                 \
        };                                                                    \
                                                                              \
        protected: inline const dictionary& settingsDict() const              \
        {                                                                     \
            return manager_.settingsDict();                                   \
        };                                                                    \
                                                                              \
        protected: virtual const dictionary& settingsRegionDict               \
        (                                                                     \
            const word& region = polyMesh::defaultRegion                      \
        ) const                                                               \
        {                                                                     \
            return manager_.settingsRegionDict(region);                       \
        };                                                                    \
                                                                              \
        protected: inline const dictionary& storageDict() const               \
        {                                                                     \
            return manager_.storageDict();                                    \
        };                                                                    \
                                                                              \
        protected: virtual const dictionary& storageRegionDict                \
        (                                                                     \
            const word& region = polyMesh::defaultRegion                      \
        ) const                                                               \
        {                                                                     \
            return manager_.storageRegionDict(region);                        \
        };                                                                    \
                                                                              \
        protected: virtual bool storageItemDict                               \
        (                                                                     \
            const word& item,                                                 \
            dictionary& dict,                                                 \
            const word& region = polyMesh::defaultRegion                      \
        ) const                                                               \
        {                                                                     \
            return manager_.storageItemDict(item, dict, region);              \
        };                                                                    \
                                                                              \
        public: NAME(const MANAGER& manager)                                  \
        :                                                                     \
            solverManager<MANAGER::Mesh>::TYPE                                \
            (                                                                 \
                manager.args(), manager.time(), manager.mesh()                \
            ),                                                                \
            manager_(manager)                                                 \
        {                                                                     \
            init();                                                           \
        };                                                                    \
                                                                              \
        public: virtual ~NAME() {};


# define solverManagerNestedClass_add(TYPE, ITEM)                             \
                                                                              \
        private: mutable autoPtr<TYPE> ITEM##Ptr_;                            \
                                                                              \
        public: TYPE& ITEM() const                                            \
        {                                                                     \
            if (ITEM##Ptr_.empty())                                           \
            {                                                                 \
                ITEM##Ptr_.set(new TYPE(this->manager()));                    \
            }                                                                 \
                                                                              \
            return ITEM##Ptr_();                                              \
        }


# define solverManagerNestedClass_addVirtual(TYPE, ITEM)                      \
                                                                              \
        private: mutable autoPtr<TYPE> ITEM##Ptr_;                            \
                                                                              \
        public: virtual TYPE& ITEM() const                                    \
        {                                                                     \
            if (ITEM##Ptr_.empty())                                           \
            {                                                                 \
                ITEM##Ptr_.set(new TYPE(this->manager()));                    \
            }                                                                 \
                                                                              \
            return ITEM##Ptr_();                                              \
        }


# define solverManagerNestedClass_end()                                       \
                                                                              \
    };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define solverManager_add(TYPE, ITEM)                                        \
                                                                              \
        private: mutable autoPtr<TYPE> ITEM##Ptr_;                            \
                                                                              \
        public: TYPE& ITEM() const                                            \
        {                                                                     \
            if (ITEM##Ptr_.empty())                                           \
            {                                                                 \
                ITEM##Ptr_.set(new TYPE(*this));                              \
            }                                                                 \
                                                                              \
            return ITEM##Ptr_();                                              \
        }


# define solverManager_addVirtual(TYPE, ITEM)                                 \
                                                                              \
        private: mutable autoPtr<TYPE> ITEM##Ptr_;                            \
                                                                              \
        public: virtual TYPE& ITEM() const                                    \
        {                                                                     \
            if (ITEM##Ptr_.empty())                                           \
            {                                                                 \
                ITEM##Ptr_.set(new TYPE(*this));                              \
            }                                                                 \
                                                                              \
            return ITEM##Ptr_();                                              \
        }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class solverManager Declaration
\*---------------------------------------------------------------------------*/

template <class MESH>
class solverManager
{
public:

    // Public typedefs

        typedef MESH Mesh;
        typedef typename MESH::Mesh MeshType;

// TODO: Macro new for regionMeshes where mesh is the mesh of one region:
//       Region and Storage should be based on MeshType
//       with a template specification for dydnamicFvMesh

        typedef solverManagerMessages<MESH> Messages;
        typedef solverManagerRegions<MESH> Regions;
        typedef solverManagerRegion<MESH> Region;
        typedef solutionControl Control;
        typedef solverManagerStorage<MESH> Storage;

    // Public Structs

// TODO: Rename
        //- Time parameters
        struct parameters
        {
            bool adjustTimeStep;
            scalar maxCo;
            scalar maxDeltaT;
            scalar CoNum;
        };


protected:

    // Protected Member Functions

        // Evolution

            //- Set time step width
            //  A return value of true means active
            virtual bool setDeltaT(scalar& deltaT) const
            {
                return false;
            };

            //- Set time parameters
            //  A return value of true means active
            virtual bool setParameters(parameters& tp) const
            {
                return false;
            };

            //- Set Courant Number
            //  A return value of true means active
            virtual bool setCoNum(scalar& CoNum) const
            {
                return false;
            };

            //- Pre loop sequence
            virtual void timePre() const {};

            //- Post loop sequence
            virtual void timePost() const {};

            //- Pre loop sequence
            virtual void loopPre() const {};

            //- Post loop sequence
            virtual void loopPost() const {};

            //- Pre run sequence
            virtual void runPre() const {};

            //- Post run sequence
            virtual void runPost() const {};

            //- Pre write sequence
            virtual void writePre() const {};

            //- Post write sequence
            virtual void writePost() const {};


private:

    // Private Data

        //- Time parameters
        mutable autoPtr<parameters> paramPtr_;

        //- Messages
        mutable autoPtr<Messages> msgPtr_;

        // Basic

            //- Reference to args
            const argList& args_;

            //- Reference to time
            Time& time_;

            //- Reference to mesh
            Mesh& mesh_;

            //- Master manager?
            const bool master_;

            //- Properties dictionary
            IOdictionary propertiesDict_;

        // Evolution

            //- Pre runTime phase?
            mutable bool prePhase_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        solverManager(const solverManager&);

        //- Disallow default bitwise assignment
        void operator=(const solverManager&);

        // Basic

            //- Error if not master manager
            void errorIfNotMaster() const;

        // Evolution

            //- Pre runTime phase?
            inline bool& prePhase() const
            {
                return prePhase_;
            };

            //- Read time parameters
            void readParameters(parameters& tc) const;

            //- Time step setting
            void calcDeltaT(scalar& deltaT) const;

            //- Time step application
            void applyDeltaT() const;


public:

    // Constructor

        solverManager
        (
            const argList& args,
            Time& time,
            Mesh& mesh,
            const word& prefix,
            const bool& master = true
        );

    // Destructor

        virtual ~solverManager() {};

    // Member Functions

        //- Time parameters
        const parameters& param() const;

        //- Messages
        const Messages& msg() const;

        //- Regions
        virtual Regions& regions() const = 0;

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

            //- Reference to mesh
            inline Mesh& mesh() const
            {
                return mesh_;
            };

            //- Is this the master manager?
            inline const bool& master() const
            {
                return master_;
            };

            //- Refrence to control dictionary
            inline const dictionary& controlDict() const
            {
                return time_.controlDict();
            }

            //- Reference to properties dictionary
            inline const dictionary& propertiesDict() const
            {
                return propertiesDict_;
            };

            //- Reference to settings dictionary
            inline const dictionary& settingsDict() const
            {
                return propertiesDict_.subDict("settings");
            };

            //- Reference to settings dictionary for named region
            const dictionary& settingsRegionDict
            (
                const word& region = polyMesh::defaultRegion
            ) const;

            //- Reference to storage dictionary
            inline const dictionary& storageDict() const
            {
                return propertiesDict_.subDict("storage");
            };

            //- Reference to storage dictionary for named region
            const dictionary& storageRegionDict
            (
                const word& region = polyMesh::defaultRegion
            ) const;

            //- Reference to dictionary of named item of storage
            bool storageItemDict
            (
                const word& item,
                dictionary& dict,
                const word& region = polyMesh::defaultRegion
            ) const;

        // Evolution

            //- Read
            virtual void read() const = 0;;

            //- Main time loop wrapper
            bool loop() const;

            //- Main time run wrapper
            bool run() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "solverManager.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
