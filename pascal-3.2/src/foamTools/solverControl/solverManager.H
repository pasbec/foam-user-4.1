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

#include "argList.H"
#include "foamTime.H"
#include "dictionary.H"
#include "solution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define createManager(NAME, MESH)                                            \
                                                                              \
    NAME##Manager manager(MESH);

# define createManagerAndData(NAME, MESH)                                     \
                                                                              \
    NAME##Manager manager(MESH);                                              \
    NAME##Manager::storage& data = manager.data();

# define setManagerScope(NAME)                                                \
                                                                              \
    NAME##Manager& NAME = manager.NAME();                                     \
    NAME##Manager& manager = NAME;

# define setManagerScopeAndData(NAME)                                         \
                                                                              \
    NAME##Manager& NAME = manager.NAME();                                     \
    NAME##Manager& manager = NAME;                                            \
    NAME##Manager::storage& data = manager.data();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define solverManagerStorage_add(TYPE, DATA)                                 \
                                                                              \
    private: TYPE* DATA##Ptr_;                                                \
    private: void create_##DATA();                                            \
    public: TYPE& DATA()                                                      \
    {                                                                         \
        if (!DATA##Ptr_)                                                      \
        {                                                                     \
            create_##DATA();                                                  \
        }                                                                     \
                                                                              \
        if (!DATA##Ptr_)                                                      \
        {                                                                     \
            FatalErrorIn(manager().type() + "::storage::" + #DATA + "()")     \
                << "Pointer for "                                             \
                    << #DATA                                                  \
                    << " is still unallocated after creation"                 \
                    << abort(FatalError);                                     \
        }                                                                     \
                                                                              \
        return *DATA##Ptr_;                                                   \
    }

# define solverManagerStorage_info(DATA)                                      \
                                                                              \
    Info << "Create "                                                         \
        << #DATA                                                              \
        << " for region " << mesh().name()                                    \
        << endl;                                                              \
                                                                              \
    if (DATA##Ptr_)                                                           \
    {                                                                         \
        FatalErrorIn(manager().type() + "::storage::create_##DATA()")         \
            << "Pointer for "                                                 \
                << #DATA                                                      \
                << " is already allocated before creation"                    \
                << abort(FatalError);                                         \
    }

# define solverManagerStorage_assert(DATA, REQUIRE)                           \
                                                                              \
    if (!REQUIRE##Ptr_)                                                       \
    {                                                                         \
        FatalErrorIn(manager().type() + "::storage::create_##DATA()")         \
            << "In order to allocate the pointer for "                        \
                << #DATA                                                      \
                << " it is necessary to allocate the pointer for "            \
                << #REQUIRE                                                   \
                << ", first"                                                  \
                << abort(FatalError);                                         \
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
protected:

    // Protected Nested Classes

        //- Messages
        class messages
        {
        private:

            // Private Data

                //- Reference to args
                const argList& args_;

                //- Reference to time
                const Time& time_;

                //- Reference to mesh
                const MESH& mesh_;

            // Private Member Functions

                //- Disallow default bitwise copy construct
                messages(const messages&);

                //- Disallow default bitwise assignment
                void operator=(const messages&);

        public:

            // Constructors

                //- Construct from mesh
                messages
                (
                    const argList& args,
                    const Time& time,
                    const MESH& mesh
                )
                :
                args_(args),
                time_(time),
                mesh_(mesh)
                {};

            // Member Functions

                //- Reference to args
                const argList& args()
                {
                    return args_;
                };

                //- Reference to time
                const Time& time()
                {
                    return time_;
                };

                //- Return reference to mesh
                const MESH& mesh() const
                {
                    return mesh_;
                }

                //- Info about starting time loop
                void newLine() const
                {
                    Info << nl;
                }

                //- Info about starting time loop
                void startTimeLoop() const
                {
                    Info << nl << "Starting time loop"
                        << nl << endl;
                }

                //- Info about current time step
                void timeStep() const
                {
                    Info << "Time step = "
                        << mesh().time().timeIndex()
                        << endl;
                };

                //- Print current time step width
                void timeDeltaT() const
                {
                    Info << "deltaT = "
                        <<  mesh().time().deltaT().value()
                        << endl;
                };

                //- Print current time
                void timeIs() const
                {
                    Info << "Time = "
                        << mesh().time().value()
                        << endl;
                };

                //- Print execution time
                void executionTime() const
                {
                    Info << "ExecutionTime = "
                        << scalar(mesh().time().elapsedCpuTime()) << " s"
                        << endl;
                };

                //- Print end statement
                void end() const
                {
                    Info << "End" << nl << endl;
                };

        };

        //- Storage
        class storage
        {
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
                storage(const storage&);

                //- Disallow default bitwise assignment
                void operator=(const storage&);

                //- Initialize data sequence
                virtual void init() = 0;

        public:

            // Constructors

                //- Construct from mesh
                storage
                (
                    const argList& args,
                    Time& time,
                    MESH& mesh
                )
                :
                args_(args),
                time_(time),
                mesh_(mesh)
                {};

            // Destructor

                virtual ~storage()
                {
                };

            // Member Functions

                //- Reference to args
                const argList& args()
                {
                    return args_;
                };

                //- Reference to time
                Time& time()
                {
                    return time_;
                };

                //- Return reference to mesh
                MESH& mesh()
                {
                    return mesh_;
                }
        };

    // Protected Data

        //- Messages
        messages* msgPtr_;

        // Basic

            //- Reference to args
            const argList& args_;

            //- Reference to time
            Time& time_;

            //- Reference to mesh
            MESH& mesh_;

            //- Master manager
            const bool master_;

            //- Properties dictionary
            IOdictionary properties_;

        // Regions

            // Default region

                //- Label of base region
                const label baseRegion_;

                //- Name of base region
                const word baseRegionName_;

                //- Solution dictionary of base region
                solution& baseSolutionDict_;

        // Evolution

    // Protected Member Functions

        //- Error if not master manager
        void errorIfNotMaster() const
        {
            if (!master_)
            {
                FatalErrorIn("solverManager::msg()")
                    << "This solverManager instantiation is NOT"
                    << " the master manager."
                    << abort(FatalError);
            }
        };

        // Evolution

            //- Courant-number calculation
            //  Needs to return true if the Courant number is beeing set
            virtual bool calcCoNum(scalar& CoNum) const
            {
                return false;
            };

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

        // Evolution

            //- Pre time-loop state
            mutable bool preRunTime_;

            //- Courant Number
            mutable scalar CoNum_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        solverManager(const solverManager&);

        //- Disallow default bitwise assignment
        void operator=(const solverManager&);

        // Evolution

            //- Refenrence to pre time-loop state
            bool& preRunTime() const
            {
                return preRunTime_;
            };

            //- Time step setting
            bool setDeltaT() const;

public:

    // Constructor

        solverManager
        (
            const argList args,
            Time& time,
            MESH& mesh,
            const word& name = "solver",
            const bool& master = true,
            const label& regionI0 = -1
        );

    // Destructor

        virtual ~solverManager()
        {
            if (master_) delete msgPtr_;
        };

    // Member Functions

        //- Messages
        const messages& msg() const
        {
            errorIfNotMaster();

            return *msgPtr_;
        };

        //- Storage
        virtual storage& data() const = 0;

        // Basic

            //- Reference to args
            const argList& args()
            {
                return args_;
            };

            //- Const reference to time
            const Time& time() const
            {
                return time_;
            };

            //- Reference to time
            Time& time()
            {
                return time_;
            };

            //- Const reference to mesh
            const MESH& mesh() const
            {
                return mesh_;
            };

            //- Reference to mesh
            MESH& mesh()
            {
                return mesh_;
            };

            //- Is this the master manager?
            const bool& master() const
            {
                return master_;
            };

            //- Reference properties dictionary
            const IOdictionary& properties() const
            {
                return properties_;
            };

        // Regions

            // Default region

                //- Reference to base region ID
                const label& base() const
                {
                    return baseRegion_;
                };

                //- Reference to base region name
                const word& baseName() const
                {
                    return baseRegionName_;
                };

                //- Reference to solution dictionary of base region
                const solution& baseSolutionDict() const
                {
                    return baseSolutionDict_;
                };

        // Evolution

            //- Return current Courant Number
            const scalar& CoNum() const
            {
                return CoNum_;
            };

            //- Main time run loop wrapper
            bool run() const;

            //- Main time write wrapper
            void write() const;

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
