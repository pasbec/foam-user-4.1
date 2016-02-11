/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

#include "foamTime.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define setControlScope(name)                                                \
                                                                              \
    name##Control& name = control.name();                                     \
    name##Control& control = name;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
public:

    // Message class
    class solverControlMessages
    {
    private:

        // Private Data

            // Basic

                //- Reference to time
                const Time& time_;

                //- Reference to region mesh
                const MESH& mesh_;

        // Private Member Functions

            //- Disallow default bitwise copy construct
            solverControlMessages(const solverControlMessages&);

            //- Disallow default bitwise assignment
            void operator=(const solverControlMessages&);



    public:

        // Constructors

            //- Construct from mesh
            solverControlMessages(const MESH& mesh)
            :
            time_(mesh.time()),
            mesh_(mesh)
            {};

        // Member Functions

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
                    << time_.timeIndex()
                    << nl << endl;
            };

            //- Print current time
            void time() const
            {
                Info << "Time = "
                    << time_.value()
                    << nl << endl;
            };

            //- Print execution time
            void executionTime() const
            {
                Info << "ExecutionTime = "
                    << scalar(time_.elapsedCpuTime()) << " s"
                    << nl << endl << endl;
            };

            //- Print end statement
            void end() const
            {
                Info << "End" << nl << endl;
            };

    };

protected:

    // Protected Data

        // Basic

            //- Reference to time
            const Time& time_;

            //- Reference to region mesh
            const MESH& mesh_;

            //- Master control
            const bool master_;

            //- Name of base region
            const word baseRegionName_;

            //- Label of base region
            const label baseRegion_;

            //- Properties dictionary
            IOdictionary propDict_;

        // Dictionaries

            //- Solution dictionary (of base region)
            IOdictionary& baseSolutionDict_;

        // Master control

            //- Message
            solverControlMessages* msg_;

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
            const MESH& mesh,
            const word& name = "solver",
            const bool& master = true,
            const label& regionI0 = -1
        );

    // Destructor

        virtual ~solverControl()
        {
            if (master_) delete msg_;
        };

    // Member Functions

        // Basic

            //- Reference to mesh
            const Time& time() const
            {
                return time_;
            };

            //- Reference to mesh
            const MESH& mesh() const
            {
                return mesh_;
            };

            //- Reference properties dictionary
            const IOdictionary& properties() const
            {
                return propDict_;
            };

        // Master control

            //- Message
            const solverControlMessages& msg() const
            {
                errorIfNotMaster();

                return *msg_;
            };

        // Properties

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

        // Time and database

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
