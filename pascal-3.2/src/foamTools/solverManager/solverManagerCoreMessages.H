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
    solverManagerCoreMessages

Description
    ...

\*---------------------------------------------------------------------------*/

#ifndef solverManagerCoreMessages_H
#define solverManagerCoreMessages_H

#include "argList.H"
#include "foamTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template <class MESH> class solverManagerCore;

/*---------------------------------------------------------------------------*\
                   Class solverManagerCoreMessages Declaration
\*---------------------------------------------------------------------------*/

template <class MESH>
class solverManagerCoreMessages
{
    // Friends
    friend class solverManagerCore<MESH>;

public:

    // Public typedefs

        typedef solverManagerCore<MESH> Manager;

        typedef MESH ManagerMesh;


protected:

    // Private data types

        //- Messages
        class Messages
        {
        private:

            // Private Data

                //- Reference to args
                const argList& args_;

                //- Reference to time
                const Time& time_;

                //- Reference to mesh
                const ManagerMesh& mesh_;

            // Private Member Functions

                //- Disallow default bitwise copy construct
                Messages(const Messages&);

                //- Disallow default bitwise assignment
                void operator=(const Messages&);


        public:

            // Constructor

                Messages
                (
                    const argList& args,
                    const Time& time,
                    const ManagerMesh& mesh
                )
                :
                    args_(args), time_(time), mesh_(mesh)
                {};

            // Destructor

                ~Messages() {};

            // Member Functions

                //- Info about starting time loop
                inline void newLine() const
                {
                    Info << nl;
                }

                //- Info about starting time loop
                inline void startTimeLoop() const
                {
                    Info << "Starting time loop"
                        << endl;
                }

                //- Info about current time step
                inline void timeStep() const
                {
                    Info << "Time step = "
                        << time_.timeIndex()
                        << endl;
                };

                //- Print current time step width
                inline void timeDeltaT() const
                {
                    Info << "deltaT = "
                        <<  time_.deltaT().value()
                        << endl;
                };

                //- Print current time
                inline void timeIs() const
                {
                    Info << "Time = "
                        << time_.value()
                        << endl;
                };

                //- Print execution time
                inline void executionTime() const
                {
                    Info << "ExecutionTime = "
                        << scalar(time_.elapsedCpuTime()) << " s"
                        << endl;
                };

                //- Print end statement
                inline void end() const
                {
                    Info << "End" << nl << endl;
                };
        };
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
