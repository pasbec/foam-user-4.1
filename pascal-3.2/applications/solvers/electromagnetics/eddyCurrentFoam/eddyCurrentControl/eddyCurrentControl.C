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

\*---------------------------------------------------------------------------*/

#include "eddyCurrentControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(eddyCurrentControl, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void eddyCurrentControl::readDictDataIfModified()
{
    // Eddy current properties
    if (propDict_.readIfModified())
    {
        // Update data
        frequency_.value() = propDict_.lookupOrDefault<dimensionedScalar>
            ("frequency", frequency_).value();
        omega_.value() = 2.0 * mathematicalConstant::pi * frequency_.value();
    }

    // Solution properties
    if
    (
        baseSolutionDict_.readIfModified()
     || conductorSolutionDict_.readIfModified()
    )
    {
        // Update subdicts
        AVdict_ = baseSolutionDict_.subDict("solvers").subDict("AV");
        Adict_ = baseSolutionDict_.subDict("solvers").subDict("A");
        Vdict_ = conductorSolutionDict_.subDict("solvers").subDict("V");

        // Update data
        tol_ = AVdict_.lookupOrDefault<scalar>("tolerance", tol_);
        tolScale_ = AVdict_.lookupOrDefault<scalar>("scale",tolScale_);
        minIter_ = AVdict_.lookupOrDefault<int>("minIter", minIter_);
        maxIter_ = AVdict_.lookupOrDefault<int>("maxIter", maxIter_);
        relax_ = AVdict_.lookupOrDefault<scalar>("relax", relax_);

        // Update subDict copies
        subAdict_ = Adict_;
        subVdict_ = Vdict_;
        subAdict_.set<scalar>("tolerance", subTol_);
        subVdict_.set<scalar>("tolerance", subTol_);
    }
}

void eddyCurrentControl::timeReset()
{
    readDictDataIfModified();

    run_ = false;
    iter_ = -1;

    subRun_ = false;
    subIter_ = -1;
    subTol_ = 1.0;
    subAdict_ = Adict_;
    subVdict_ = Vdict_;
}

void eddyCurrentControl::decreaseSubTolerance() const
{
    if ((tol_ <= subTol_) && !firstIteration())
    {
        // Scale sub-tolerance
        subTol_ *= min(mag(tolScale_), 1.0);

        // If residual is already smaller
        // adjust sub-tolerance accordingly
        subTol_ = min(subTol_, AVres_);

        // Make sure sub-tolerance does
        // not get smaller then tolerance
        subTol_ = max(subTol_, tol_);
    }

    Info << nl << "Tolerance = " << subTol_
        << " / " << tol_ << endl << nl;

    subAdict_.set<scalar>("tolerance", subTol_);
    subVdict_.set<scalar>("tolerance", subTol_);

    subIter_ = 0;
}

const bool& eddyCurrentControl::subRun() const
{
    subRun_ = !checkSubConvergence();

    subIter_++;

    return subRun_;
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

eddyCurrentControl::eddyCurrentControl(const regionControl& rfvc)
:
    regIOobject
    (
        IOobject
        (
            "eddyCurrentControl",
            rfvc.rmesh().time().timeName(),
            "uniform",
            rfvc.rmesh().time().db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    ),
    control_(rfvc),
    rmesh_(rfvc.rmesh()),
    time_(rmesh_.time()),
    propDict_
    (
        IOdictionary
        (
            IOobject
            (
                "eddyCurrentProperties",
                time_.constant(),
                time_.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    ),
    baseRegionName_
    (
        word(propDict_.subDict("regions").lookup("base"))
    ),
    conductorRegionName_
    (
        word(propDict_.subDict("regions").lookup("conductor"))
    ),
    baseRegion_
    (
        rmesh_.regionIndex(baseRegionName_)
    ),
    conductorRegion_
    (
        rmesh_.regionIndex(conductorRegionName_)
    ),
    frequency_
    (
        "frequency",
        dimensionedScalar(propDict_.lookup("frequency"))
    ),
    omega_
    (
        "omega",
        2.0 * mathematicalConstant::pi * frequency_
    ),
    baseSolutionDict_(rmesh_[baseRegion_].solutionDict()),
    conductorSolutionDict_(rmesh_[conductorRegion_].solutionDict()),
    AVdict_(baseSolutionDict_.subDict("solvers").subDict("AV")),
    Adict_(baseSolutionDict_.subDict("solvers").subDict("A")),
    Vdict_(conductorSolutionDict_.subDict("solvers").subDict("V")),
    tol_(AVdict_.lookupOrDefault<scalar>("tolerance", 1e-04)),
    tolScale_(AVdict_.lookupOrDefault<scalar>("scale", 0.1)),
    minIter_(AVdict_.lookupOrDefault<int>("minIter", 0)),
    maxIter_(AVdict_.lookupOrDefault<int>("maxIter", 100)),
    relax_(AVdict_.lookupOrDefault<scalar>("relax", 1.0)),
    stop_(false),
    run_(false),
    iter_(-1),
    subRun_(false),
    subIter_(-1),
    subTol_(1.0),
    subAdict_(Adict_),
    subVdict_(Vdict_),
    AVres_(subTol_-SMALL),
    Ares_(AVres_),
    Vres_(AVres_)
{
    // TODO: Check tol_ < 1 ???
    // TODO: Check maxIter_ > 0
    // TODO: Check maxIter_ > minIter_
    // TODO: Check 0 < relax_ <= 1
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void eddyCurrentControl::relax(volVectorField& vf) const
{
    if (checkRelax())
    {
        if (!firstIteration())
        {
            Info << "Relax " << vf.name()
                << ", Factor = " << relax_ << endl;

            vf.relax(relax_);
        }

        vf.storePrevIter();
    }
}

void eddyCurrentControl::relax(const regionVolVectorField& rvf) const
{
    if (checkRelax())
    {
        if (!firstIteration())
        {
            Info << "Relax " << rvf.name()
                << " = " << relax_ << endl;

            forAll(rvf.mesh().regionNames(), regionI)
            {
                rvf[regionI].relax(relax_);
            }
        }

        forAll(rvf.mesh().regionNames(), regionI)
        {
            rvf[regionI].storePrevIter();
        }
    }
}

const bool& eddyCurrentControl::run()
{
    bool converged = checkConvergence();
    bool subConverged = checkSubConvergence();
//     bool iterAboveMin = checkMinIterations(); // TODO: Makes only sense with relTol
    bool iterBelowMax = checkMaxIterations();

    // Update data if last subiteration has converged
    if (!converged && subConverged)
    {
        readDictDataIfModified();
    }

    // Run if not converged and below max iterations
    run_ = !converged && iterBelowMax;

    if (run_)
    {
        iter_++;

        if (!subRun())
        {
            decreaseSubTolerance();
        }
    }
    else
    {
        timeReset();
    }

    return run_;
}


void eddyCurrentControl::setResidualOfA(scalar residual) const
{
    Ares_ = residual;

    updateResidual();
}

void eddyCurrentControl::setResidualOfV(scalar residual) const
{
    Vres_ = residual;

    updateResidual();
}


void eddyCurrentControl::subWrite() const
{
    if
    (
        checkSubConvergence()
     && !checkConvergence()
    )
    {
        Info << "Write current solution" << endl;

        Time& time = const_cast<Time&>(time_);

        // TODO: Why is writeNow() not const,
        //       but write() is?
        time.writeNow();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

