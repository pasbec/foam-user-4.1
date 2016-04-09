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

template <class MESH, unsigned int NREGIONS>
void solverManager<MESH, NREGIONS>::errorIfNotMaster() const
{
    if (!master())
    {
        FatalErrorIn("solverManager::messages()")
            << "This solverManager instantiation is NOT"
            << " the master manager."
            << abort(FatalError);
    }
}


template <class MESH, unsigned int NREGIONS>
void solverManager<MESH, NREGIONS>::readParameters
(
    Parameters& parameters
) const
{
    errorIfNotMaster();

    readTimeParameters(parameters.time);
}


template <class MESH, unsigned int NREGIONS>
void solverManager<MESH, NREGIONS>::readTimeParameters
(
    TimeParameters& timePar
) const
{
    errorIfNotMaster();

    timePar.adjustTimeStep =
        controlDict().lookupOrDefault("adjustTimeStep", false);

    timePar.maxCo =
        controlDict().lookupOrDefault("maxCo", scalar(1.0));

    timePar.maxDeltaT =
        controlDict().lookupOrDefault("maxDeltaT", scalar(GREAT));

    timePar.CoNum = 0.0;
}


template <class MESH, unsigned int NREGIONS>
void solverManager<MESH, NREGIONS>::calcDeltaT
(
    scalar& deltaT
) const
{
    errorIfNotMaster();

    Parameters& par = const_cast<Parameters&>(parameters());
    TimeParameters& timePar = par.time;

    Time& runTime = const_cast<Time&>(time());

    if (timePar.adjustTimeStep)
    {

        scalar maxDeltaTFact = timePar.maxCo/(timePar.CoNum + SMALL);
        scalar deltaTFact =
            min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

        deltaT =
            min
            (
                deltaTFact * runTime.deltaT().value(),
                timePar.maxDeltaT
            );
    }
}


template <class MESH, unsigned int NREGIONS>
void solverManager<MESH, NREGIONS>::applyDeltaT() const
{
    errorIfNotMaster();

    Parameters& par = const_cast<Parameters&>(parameters());
    TimeParameters& timePar = par.time;

    Time& runTime = const_cast<Time&>(time());

    scalar deltaT = GREAT;

    bool deltaTset = setDeltaT(deltaT);

    if (!deltaTset)
    {
        bool parametersSet = setTimeParameters(timePar);

        if (!parametersSet)
        {
            readTimeParameters(timePar);
        }

        bool CoNumSet = setCoNum(timePar.CoNum);

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

        messages().newLine();

        calcDeltaT(deltaT);
    }

    runTime.setDeltaT(deltaT);

    messages().timeDeltaT();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class MESH, unsigned int NREGIONS>
solverManager<MESH, NREGIONS>::solverManager
(
    const argList& args,
    Time& time,
    ManagerMesh& mesh,
    const word& prefix,
    const bool& master
)
:
    solverManagerName(),
    parametersPtr_(NULL),
    messagesPtr_(NULL),
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

template <class MESH, unsigned int NREGIONS>
const typename solverManager<MESH, NREGIONS>::Parameters&
solverManager<MESH, NREGIONS>::parameters() const
{
    errorIfNotMaster();

    if (parametersPtr_.empty())
    {
        parametersPtr_.set(new Parameters());
        {
            readParameters(parametersPtr_());
        }
    }

    return parametersPtr_();
}


template <class MESH, unsigned int NREGIONS>
const typename solverManager<MESH, NREGIONS>::Messages&
solverManager<MESH, NREGIONS>::messages() const
{
    errorIfNotMaster();

    if (messagesPtr_.empty())
    {
        messagesPtr_.set
        (
            new Messages(args(), time(), mesh())
        );
    }

    return messagesPtr_();
}


template <class MESH, unsigned int NREGIONS>
void solverManager<MESH, NREGIONS>::read() const
{
    settings().checkRead();
    regions().checkRead();
}


template <class MESH, unsigned int NREGIONS>
void solverManager<MESH, NREGIONS>::init() const
{
    storage().checkInit();
    regions().checkInit();
}


template <class MESH, unsigned int NREGIONS>
bool solverManager<MESH, NREGIONS>::once() const
{
    errorIfNotMaster();

    Time& runTime = const_cast<Time&>(time());

    if (prePhase())
    {
        init();

        if (!args().optionFound("overwrite"))
        {
            runTime++;
        }

        messages().timeIs();
        messages().newLine();

        prePhase() = false;

        return true;
    }
    else
    {
        next();

        time().write();

        if (time().outputTime())
        {
            write();
        }

        messages().executionTime();
        messages().newLine();

        finalize();

        messages().end();

        return false;
    }
}


template <class MESH, unsigned int NREGIONS>
bool solverManager<MESH, NREGIONS>::loop() const
{
    errorIfNotMaster();

    Time& runTime = const_cast<Time&>(time());

    if (prePhase())
    {
        init();

        messages().newLine();
        messages().startTimeLoop();
        messages().newLine();

        prePhase() = false;
    }
    else
    {
        next();

        time().write();

        if (time().outputTime())
        {
            write();
        }

        messages().executionTime();
        messages().newLine();
    }

    if (runTime.loop())
    {
        read();

        messages().timeIs();
        messages().newLine();

        return true;
    }
    else
    {
        finalize();

        messages().end();

        return false;
    }
}


template <class MESH, unsigned int NREGIONS>
bool solverManager<MESH, NREGIONS>::run() const
{
    errorIfNotMaster();

    Time& runTime = const_cast<Time&>(time());

    if (prePhase())
    {
        init();

        messages().newLine();
        messages().startTimeLoop();
        messages().newLine();

        prePhase() = false;
    }
    else
    {
        next();

        time().write();

        if (time().outputTime())
        {
            write();
        }

        messages().executionTime();
        messages().newLine();
    }

    if (runTime.run())
    {
        read();

        applyDeltaT();

        runTime++;

        messages().timeIs();
        messages().newLine();

        return true;
    }
    else
    {
        finalize();

        messages().end();

        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

