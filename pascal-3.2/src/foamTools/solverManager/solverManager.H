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
                                                                              \
        private: const MANAGER& manager_;                                     \
                                                                              \
        private: NAME(const NAME&);                                           \
                                                                              \
        private: void operator=(const NAME&);                                 \
                                                                              \
        protected: inline const MANAGER& manager() const                      \
        {                                                                     \
            return manager_;                                                  \
        };                                                                    \
                                                                              \
        protected: inline const dictionary& properties() const                \
        {                                                                     \
            return manager_.properties();                                     \
        };                                                                    \
                                                                              \
        protected: inline const dictionary& modules() const                   \
        {                                                                     \
            return manager_.modules();                                        \
        };                                                                    \
                                                                              \
        protected: inline bool modules                                        \
        (                                                                     \
            const word& module,                                               \
            dictionary& dict,                                                 \
            const Switch& negate = false                                      \
        ) const                                                               \
        {                                                                     \
            return manager_.modules(module, dict, negate);                    \
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
            MESH& mesh_;

            //- Master manager?
            const bool master_;

            //- Properties dictionary
            IOdictionary properties_;

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
            MESH& mesh,
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

        // Regions

            //- All regions
            virtual Regions& regions() const = 0;

            //- Default region
            virtual Region& defaultRegion() const = 0;

        //- Control
        virtual Control& control() const = 0;

        //- Storage
        virtual Storage& storage() const = 0;


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

            //- Refrence to controlDict
            inline const dictionary& controlDict() const
            {
                return time_.controlDict();
            }

            //- Reference to mesh
            inline MESH& mesh() const
            {
                return mesh_;
            };

            //- Is this the master manager?
            inline const bool& master() const
            {
                return master_;
            };

            //- Reference properties dictionary
            inline const dictionary& properties() const
            {
                return properties_;
            };

            //- Reference modules sub-dictionary
            inline const dictionary& modules() const
            {
                return properties_.subDict("modules");
            };

// TODO: Move to C-file
            //- Reference to named module sub-dictionary
            bool modules
            (
                const word& module,
                dictionary& dict,
                const Switch& negate = false
            ) const
            {
                if
                (
                    modules().found(module)
                 && modules().isDict(module)
                )
                {
                    dict = modules().subDict(module);

                    Switch enabled;
                    if (dict.readIfPresent("enabled", enabled))
                    {
                        return enabled;
                    }
                }

                if (!negate)
                {
                    return false;
                }
                else
                {
                    return true;
                }
            };
        // Evolution

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
