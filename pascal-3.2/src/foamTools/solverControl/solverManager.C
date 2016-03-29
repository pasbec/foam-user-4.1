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

\*---------------------------------------------------------------------------*/

#include "solverManager.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class MESH>
bool solverManager<MESH>::setDeltaT() const
{
    errorIfNotMaster();

    Time& runTime = const_cast<Time&>(time());

    const dictionary& controlDict = runTime.controlDict();

    Switch adjustTimeStep
    (
        controlDict.lookup("adjustTimeStep")
    );

    if (adjustTimeStep)
    {
        scalar maxCo
        (
            readScalar(controlDict.lookup("maxCo"))
        );

        scalar maxDeltaT =
            controlDict.lookupOrDefault<scalar>("maxDeltaT", GREAT);

        scalar maxDeltaTFact = maxCo/(CoNum() + SMALL);
        scalar deltaTFact =
            min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

        runTime.setDeltaT
        (
            min
            (
                deltaTFact * runTime.deltaT().value(),
                maxDeltaT
            )
        );

        msg().timeDeltaT();

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class MESH>
solverManager<MESH>::solverManager
(
    const argList args,
    Time& time,
    MESH& mesh,
    const word& name,
    const bool& master,
    const label& regionI0
)
:
    msgPtr_(NULL),
    args_(args),
    time_(time),
    mesh_(mesh),
    master_(master),
    properties_
    (
        IOdictionary
        (
            IOobject
            (
                name + "Properties",
                time_.constant(),
                time_.db(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        )
    ),
    baseRegion_(regionI0),
    baseRegionName_(polyMesh::defaultRegion),
    baseSolutionDict_
    (
        const_cast<solution&>
        (
            dynamicCast<const solution&>(mesh_.solutionDict())
        )
    ),
    preRunTime_(true),
    CoNum_(0.0)
{
    if (master_) msgPtr_ = new messages(args_, time_, mesh_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class MESH>
bool solverManager<MESH>::run() const
{
    errorIfNotMaster();

    Time& runTime = const_cast<Time&>(time());

    if (preRunTime())
    {
        msg().startTimeLoop();

        preRunTime() = false;
    }

    if (runTime.run())
    {
        runPre();

        bool CoNumSet = calcCoNum(CoNum_);

        if (!CoNumSet)
        {
            FatalErrorIn("solverManager::run()")
                << "A Courant Number needs to be calculated for run(). "
                << "Use the virtual function calcCoNum(scalar& CoNum) to "
                << "calulate it, which finally needs to return true."
                << abort(FatalError);
        }

        setDeltaT();

        runTime++;

        msg().timeIs();

        runPost();

        return true;
    }
    else
    {
        msg().end();

        return false;
    }
}


template <class MESH>
void solverManager<MESH>::write() const
{
    errorIfNotMaster();

    {
        writePre();

        time().write();

        writePost();

        msg().executionTime();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

