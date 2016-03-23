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
    solverControl

Description
    ...

SourceFiles
    solverControl.C

\*---------------------------------------------------------------------------*/

#ifndef solverControl_H
#define solverControl_H

#include "argList.H"
#include "foamTime.H"
#include "IOdictionary.H"
#include "solution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define createControl(NAME, MESH)                                            \
                                                                              \
    NAME##Control control(MESH);

# define createControlAndData(NAME, MESH)                                     \
                                                                              \
    NAME##Control control(MESH);                                              \
    NAME##Control::storage& data = control.data();

# define setControlScope(NAME)                                                \
                                                                              \
    NAME##Control& NAME = control.NAME();                                     \
    NAME##Control& control = NAME;

# define setControlScopeWithData(NAME)                                        \
                                                                              \
    NAME##Control& NAME = control.NAME();                                     \
    NAME##Control& control = NAME;                                            \
    NAME##Control::storage& data = control.data();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define solverControlStorage_add(TYPE, DATA)                                 \
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
            FatalErrorIn(control().type() + "::storage::" + #DATA + "()")     \
                << "Pointer for "                                             \
                    << #DATA                                                  \
                    << " is still unallocated after creation"                 \
                    << abort(FatalError);                                     \
        }                                                                     \
                                                                              \
        return *DATA##Ptr_;                                                   \
    }

# define solverControlStorage_info(DATA)                                      \
                                                                              \
    Info << "Create "                                                         \
        << #DATA                                                              \
        << " for region " << mesh().name()                                    \
        << endl;                                                              \
                                                                              \
    if (DATA##Ptr_)                                                           \
    {                                                                         \
        FatalErrorIn(control().type() + "::storage::create_##DATA()")         \
            << "Pointer for "                                                 \
                << #DATA                                                      \
                << " is already allocated before creation"                    \
                << abort(FatalError);                                         \
    }

# define solverControlStorage_assert(DATA, REQUIRE)                           \
                                                                              \
    if (!REQUIRE##Ptr_)                                                       \
    {                                                                         \
        FatalErrorIn(control().type() + "::storage::create_##DATA()")         \
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
                       Class solverControl Declaration
\*---------------------------------------------------------------------------*/

template <class MESH>
class solverControl
:
    public IOdictionary
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
                        << nl << endl;
                };

                //- Print current time
                void timeIs() const
                {
                    Info << "Time = "
                        << mesh().time().value()
                        << nl << endl;
                };

                //- Print execution time
                void executionTime() const
                {
                    Info << "ExecutionTime = "
                        << scalar(mesh().time().elapsedCpuTime()) << " s"
                        << nl << endl << endl;
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
                virtual void init() {};

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

        //- Basic

            //- Reference to args
            const argList& args_;

            //- Reference to time
            Time& time_;

            //- Reference to mesh
            MESH& mesh_;

            //- Master control
            const bool master_;

            //- Name of base region
            const word baseRegionName_;

            //- Label of base region
            const label baseRegion_;

            //- Properties dictionary
            IOdictionary propDict_;

        //- Dictionaries

            //- Solution dictionary of base region
            solution& baseSolutionDict_;

    // Protected Member Functions

        //- Error if not master control
        void errorIfNotMaster() const
        {
            if (!master_)
            {
                FatalErrorIn("solverControl::msg()")
                    << "This solverControl instantiation is NOT"
                    << " the master control."
                    << abort(FatalError);
            }
        };

private:

    // Private Data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        solverControl(const solverControl&);

        //- Disallow default bitwise assignment
        void operator=(const solverControl&);

public:

    // Constructors

        //- Construct from mesh and name
        solverControl
        (
            const argList args,
            Time& time,
            MESH& mesh,
            const word& name = "solver",
            const bool& master = true,
            const label& regionI0 = -1
        );

    // Destructor

        virtual ~solverControl()
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
        virtual storage& data() = 0;

        //- Basic

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

            //- Reference to mesh
            MESH& mesh()
            {
                return mesh_;
            };

            //- Is this master control?
            const bool& master() const
            {
                return master_;
            };

            //- Reference properties dictionary
            const IOdictionary& properties() const
            {
                return propDict_;
            };

        //- Properties

            //- Reference to base region name
            const word& baseName() const
            {
                return baseRegionName_;
            };

            //- Reference to base region ID
            const label& base() const
            {
                return baseRegion_;
            };

        //- Time and database

            //- WriteData member function required by regIOobject
            bool writeData(Ostream&) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "solverControl.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
