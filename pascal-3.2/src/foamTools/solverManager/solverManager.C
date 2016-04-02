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
void solverManager<MESH>::errorIfNotMaster() const
{
    if (!master())
    {
        FatalErrorIn("solverManager::msg()")
            << "This solverManager instantiation is NOT"
            << " the master manager."
            << abort(FatalError);
    }
}


template <class MESH>
void solverManager<MESH>::readParameters
(
    parameters& tp
) const
{
    errorIfNotMaster();

    tp.adjustTimeStep =
        controlDict().lookupOrDefault("adjustTimeStep", false);

    tp.maxCo =
        controlDict().lookupOrDefault("maxCo", scalar(1.0));

    tp.maxDeltaT =
        controlDict().lookupOrDefault("maxDeltaT", scalar(GREAT));
}


template <class MESH>
void solverManager<MESH>::calcDeltaT
(
    scalar& deltaT
) const
{
    errorIfNotMaster();

    parameters& tp = const_cast<parameters&>(param());

    Time& runTime = const_cast<Time&>(time());

    if (tp.adjustTimeStep)
    {

        scalar maxDeltaTFact = tp.maxCo/(tp.CoNum + SMALL);
        scalar deltaTFact =
            min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

        deltaT =
            min
            (
                deltaTFact * runTime.deltaT().value(),
                tp.maxDeltaT
            );
    }
}


template <class MESH>
void solverManager<MESH>::applyDeltaT() const
{
    errorIfNotMaster();

    parameters& tp = const_cast<parameters&>(param());

    Time& runTime = const_cast<Time&>(time());

    scalar deltaT = GREAT;

    bool deltaTset = setDeltaT(deltaT);

    if (!deltaTset)
    {
        bool parametersSet = setParameters(tp);

        if (!parametersSet)
        {
            readParameters(tp);
        }

        bool CoNumSet = setCoNum(tp.CoNum);

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

        msg().newLine();

        calcDeltaT(deltaT);
    }

    runTime.setDeltaT(deltaT);

    msg().timeDeltaT();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class MESH>
solverManager<MESH>::solverManager
(
    const argList& args,
    Time& time,
    ManagerMesh& mesh,
    const word& prefix,
    const bool& master
)
:
    paramPtr_(NULL),
    msgPtr_(NULL),
    args_(args),
    time_(time),
    mesh_(mesh),
    master_(master),
    propertiesDict_
    (
        IOdictionary
        (
            IOobject
            (
                prefix + "Properties",
                time_.constant(),
                time_.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    ),
    prePhase_(true)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class MESH>
const typename solverManager<MESH>::parameters& solverManager<MESH>::param() const
{
    errorIfNotMaster();

    if (paramPtr_.empty())
    {
        paramPtr_.set(new parameters());
        {
            readParameters(paramPtr_());
            paramPtr_->CoNum = 0.0;
        }
    }

    return paramPtr_();
}


template <class MESH>
const typename solverManager<MESH>::Messages& solverManager<MESH>::msg() const
{
    errorIfNotMaster();

    if (msgPtr_.empty())
    {
        msgPtr_.set
        (
            new Messages(args(), time(), mesh())
        );
    }

    return msgPtr_();
}


template <class MESH>
const dictionary& solverManager<MESH>::settingsRegionDict
(
    const word& region
) const
{
    return settingsDict().subDict(region);
};


template <class MESH>
const dictionary& solverManager<MESH>::storageRegionDict
(
    const word& region
) const
{
    return storageDict().subDict(region);
};


template <class MESH>
bool solverManager<MESH>::storageItemDict
(
    const word& item,
    dictionary& dict,
    const word& region
) const
{
    const dictionary& regionDict = storageRegionDict(region);

    if
    (
        regionDict.found(item)
     && regionDict.isDict(item)
    )
    {
        dict = regionDict.subDict(item);

        Switch enabled;

        if (dict.readIfPresent("enabled", enabled))
        {
            return enabled;
        }
    }

    return false;
}


template <class MESH>
bool solverManager<MESH>::loop() const
{
    errorIfNotMaster();

    Time& runTime = const_cast<Time&>(time());

    if (prePhase())
    {
        init();

        msg().newLine();
        msg().startTimeLoop();
        msg().newLine();

        prePhase() = false;
    }
    else
    {
        writePre();

        time().write();

        writePost();

        msg().executionTime();
        msg().newLine();
    }

    if (runTime.loop())
    {
        read();

        msg().timeIs();
        msg().newLine();

        return true;
    }
    else
    {
        finalize();

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
        init();

        msg().newLine();
        msg().startTimeLoop();
        msg().newLine();

        prePhase() = false;
    }
    else
    {
        writePre();

        time().write();

        writePost();

        msg().executionTime();
        msg().newLine();
    }

    if (runTime.run())
    {
        read();

        applyDeltaT();

        runTime++;

        msg().timeIs();
        msg().newLine();

        return true;
    }
    else
    {
        finalize();

        msg().end();

        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

