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
#include "dimensionSet.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class MESH>
void solverManager<MESH>::readControls
(
    timeControls& tc
) const
{
    errorIfNotMaster();

    tc.adjustTimeStep =
        controlDict().lookupOrDefault("adjustTimeStep", false);

    tc.maxCo =
        controlDict().lookupOrDefault("maxCo", scalar(1.0));

    tc.maxDeltaT =
        controlDict().lookupOrDefault("maxDeltaT", scalar(GREAT));
}


template <class MESH>
void solverManager<MESH>::calcDeltaT
(
    scalar& deltaT
) const
{
    errorIfNotMaster();

    timeControls& tc = const_cast<timeControls&>(controls());

    Time& runTime = const_cast<Time&>(time());

    if (tc.adjustTimeStep)
    {

        scalar maxDeltaTFact = tc.maxCo/(tc.CoNum + SMALL);
        scalar deltaTFact =
            min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

        deltaT =
            min
            (
                deltaTFact * runTime.deltaT().value(),
                tc.maxDeltaT
            );
    }
}


template <class MESH>
void solverManager<MESH>::applyDeltaT() const
{
    errorIfNotMaster();

    timeControls& tc = const_cast<timeControls&>(controls());

    Time& runTime = const_cast<Time&>(time());

    scalar deltaT = GREAT;

    bool deltaTset = setDeltaT(deltaT);

    if (!deltaTset)
    {
        bool controlsSet = setControls(tc);
    
        if (!controlsSet)
        {
            readControls(tc);
        }
    
        bool CoNumSet = setCoNum(tc.CoNum);
    
        if (!CoNumSet)
        {
            FatalErrorIn("solverManager::applyDeltaT()")
                << "A Courant Number needs to be calculated if "
                    << "solverManager::setDeltaT(scalar& deltaT) "
                    << "is not beeing used. Utilize the virtual function "
                    << "setCoNum(scalar& newCoNum) to set it. The "
                    << "return value needs to be true."
                    << abort(FatalError);
        }
    
        calcDeltaT(deltaT);
    }

    runTime.setDeltaT(deltaT);

    msg().timeDeltaT();
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
    timeControlsPtr_(NULL),
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
    prePhase_(true)
{
    if (master_)
    {
        timeControlsPtr_ = new timeControls();
        {
            readControls(*timeControlsPtr_);
            timeControlsPtr_->CoNum = 0.0;
        }

        msgPtr_ = new messages(args_, time_, mesh_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class MESH>
bool solverManager<MESH>::loop() const
{
    errorIfNotMaster();

    Time& runTime = const_cast<Time&>(time());

    if (prePhase())
    {
        timePre();

        msg().startTimeLoop();

        prePhase() = false;
    }

    loopPre();

    if (runTime.loop())
    {
        msg().timeIs();

        loopPost();

        return true;
    }
    else
    {
        timePost();

        msg().end();

        return false;
    }
}

template <class MESH>
bool solverManager<MESH>::run() const
{
    errorIfNotMaster();

    Time& runTime = const_cast<Time&>(time());

    if (prePhase())
    {
        timePre();

        msg().startTimeLoop();

        prePhase() = false;
    }

    runPre();

    if (runTime.run())
    {
        applyDeltaT();

        runTime++;

        msg().timeIs();

        runPost();

        return true;
    }
    else
    {
        timePost();

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

