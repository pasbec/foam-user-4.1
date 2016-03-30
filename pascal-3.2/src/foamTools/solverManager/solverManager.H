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
#include "solutionControl.H"

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

# define solverManagerStorage_new(MANAGER)                                    \
                                                                              \
    private: storage(const storage&);                                         \
                                                                              \
    private: void operator=(const storage&);                                  \
                                                                              \
    public: storage(const MANAGER& manager)                                   \
    :                                                                         \
        solverManager<dynamicFvMesh>::storage                                 \
        (                                                                     \
            manager.args(), manager.time(), manager.mesh(),                   \
            manager.typeName + word("Storage")                                \
        )                                                                     \
    {                                                                         \
        init();                                                               \
    };                                                                        \
                                                                              \
    public: virtual ~storage() {};                                            \
                                                                              \
    public: virtual void init() const;

# define solverManagerStorage_add(TYPE, DATA)                                 \
                                                                              \
    private: mutable autoPtr<TYPE> DATA##Ptr_;                                \
                                                                              \
    private: void create_##DATA(const word init = "default") const;           \
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
        if (DATA##IsEmpty())                                                  \
        {                                                                     \
            create_##DATA(init);                                              \
        }                                                                     \
                                                                              \
        if (DATA##IsEmpty())                                                  \
        {                                                                     \
            FatalErrorIn(typeName + #DATA + "()")                             \
                << "Pointer for "                                             \
                    << #DATA                                                  \
                    << " is still unallocated after creation"                 \
                    << abort(FatalError);                                     \
        }                                                                     \
                                                                              \
        return DATA##Ptr_();                                                  \
    };                                                                        \
                                                                              \
    public: TYPE* DATA##Ptr(const word init = "default") const                \
    {                                                                         \
        if (DATA##IsEmpty())                                                  \
        {                                                                     \
            create_##DATA(init);                                              \
        }                                                                     \
                                                                              \
        if (DATA##IsEmpty())                                                  \
        {                                                                     \
            FatalErrorIn(typeName + #DATA + "()")                             \
                << "Pointer for "                                             \
                    << #DATA                                                  \
                    << " is still unallocated after creation"                 \
                    << abort(FatalError);                                     \
        }                                                                     \
                                                                              \
        return DATA##Ptr_.operator->();                                       \
    };

# define solverManagerStorage_info(DATA)                                      \
                                                                              \
    Info << "Create "                                                         \
        << #DATA                                                              \
        << " for region " << mesh().name()                                    \
        << endl;                                                              \
                                                                              \
    if (!DATA##IsEmpty())                                                     \
    {                                                                         \
        FatalErrorIn(typeName + "::create_" + #DATA + "()")                   \
            << "Pointer for "                                                 \
                << #DATA                                                      \
                << " is already allocated before creation"                    \
                << abort(FatalError);                                         \
    }

# define solverManagerStorage_assert(DATA, REQUIRE)                           \
                                                                              \
    if (!REQUIRE##IsValid())                                                  \
    {                                                                         \
        FatalErrorIn(typeName + "::create_" + #DATA + "()")                   \
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
public:

    // Public structs

//         //- Time parameters
        struct parameters
        {
            bool adjustTimeStep;
            scalar maxCo;
            scalar maxDeltaT;
            scalar CoNum;
        };


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

            // Constructor

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

            // Destructor

                ~messages() {};

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
                    Info << "Starting time loop"
                        << endl;
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


        public:

            const word typeName;

            // Constructor

                //- Construct from mesh
                storage
                (
                    const argList& args,
                    Time& time,
                    MESH& mesh,
                    const word& name
                )
                :
                args_(args),
                time_(time),
                mesh_(mesh),
                typeName(name)
                {};

            // Destructor

                virtual ~storage() {};

            // Member Functions

                //- Reference to args
                const argList& args() const
                {
                    return args_;
                };

                //- Reference to time
                Time& time() const
                {
                    return time_;
                };

                //- Return reference to mesh
                MESH& mesh() const
                {
                    return mesh_;
                }

                //- Initialize data
                virtual void init() const = 0;
        };

    // Protected Data

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
        mutable parameters* paramPtr_;

        //- Messages
        mutable messages* msgPtr_;

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

            //- Pre runTime phase?
            mutable bool prePhase_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        solverManager(const solverManager&);

        //- Disallow default bitwise assignment
        void operator=(const solverManager&);

        //- Pre runTime phase
        bool& prePhase() const
        {
            return prePhase_;
        };

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
            if (master_)
            {
                delete paramPtr_;
                delete msgPtr_;
            }
        };

    // Member Functions

        //- Time parameters
        const parameters& param() const
        {
            errorIfNotMaster();

            return *paramPtr_;
        }

        //- Messages
        const messages& msg() const
        {
            errorIfNotMaster();

            return *msgPtr_;
        };

        //- Control
        virtual solutionControl& control() const = 0;

        //- Storage
        virtual storage& data() const = 0;

        // Basic

            //- Reference to args
            const argList& args() const
            {
                return args_;
            };

            //- Reference to time
            Time& time() const
            {
                return time_;
            };

            //- Refrence to controlDict
            const dictionary& controlDict() const
            {
                return time_.controlDict();
            }

            //- Reference to mesh
            MESH& mesh() const
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

            //- Initialize solver manager
            void init() const;

            //- Main time loop wrapper
            bool loop() const;

            //- Main time run wrapper
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
