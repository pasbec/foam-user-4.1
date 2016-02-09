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
        relTol_ = AVdict_.lookupOrDefault<scalar>("relTol", relTol_);
        outerRelTol_ = AVdict_.lookupOrDefault<int>("outerRelTol", outerRelTol_);
        outerMaxIter_ = AVdict_.lookupOrDefault<int>("outerMaxIter", outerMaxIter_);
        relax_ = AVdict_.lookupOrDefault<scalar>("relax", relax_);

        // Reset sub-dictionaries
        resetSubDictionaries();
    }
}

//- Check convergence state
void eddyCurrentControl::updateResidual() const
{
    if (!oldAVresSet_)
    {
        oldAVres_ = AVres_;
        oldAVrelRes_ = AVrelRes_;
        oldAVresSet_ = true;
    }

    AVres_ = returnReduce(max(Ares_, Vres_), maxOp<scalar>());
    AVrelRes_ = mag(oldAVres_ - AVres_)/(oldAVres_ + VSMALL);
};

void eddyCurrentControl::resetSubDictionaries() const
{
    // Create copy
    subAdict_ = Adict_;
    subVdict_ = Vdict_;

    // Override tolerance
    subAdict_.set<scalar>
    (
        "tolerance",
        min(tol_, Adict_.lookupOrDefault<scalar>("tolerance", tol_))
    );
    subVdict_.set<scalar>
    (
        "tolerance",
        min(tol_, Vdict_.lookupOrDefault<scalar>("tolerance", tol_))
    );

    // Override relative tolerance
    subAdict_.set<scalar>("relTol", relTol_);
    subVdict_.set<scalar>("relTol", relTol_);
}

void eddyCurrentControl::decreaseSubTolerance() const
{
    // Linear convergence progress
    scalar x = min(max(AVres_-tol_, 0.0), 1.0);

    // Relative tolerance decrease function (atan, scaled)
    scalar decRelTol = atan(subScale_ * x) / mathematicalConstant::piByTwo;

    // Relax decrease function (atan)
    scalar decRelax = atan(10*x) / mathematicalConstant::piByTwo;

    subRelTol_ = relTol_ * decRelTol;
    subAdict_.set<scalar>("relTol", subRelTol_);
    subVdict_.set<scalar>("relTol", subRelTol_);

    subRelax_ = relax_ + (1 - relax_) * (1 - decRelax);

//     Info << endl;
//     Info << "loop_ = " << loop_ << endl;
//     Info << "iter_ = " << iter_ << endl;
//     Info << "AVres_ = " << AVres_ << endl;
//     Info << "AVrelRes_ = " << AVrelRes_ << endl;
//     Info << "oldAVres_ = " << oldAVres_ << endl;
//     Info << "oldAVrelRes_ = " << oldAVrelRes_ << endl;
//
//     Info << "subLoop_ = " << subLoop_ << endl;
//     Info << "subIter_ = " << subIter_ << endl;
//     Info << "subRelTol_ = " << subRelTol_ << endl;
//     Info << "subScale_ = " << subScale_ << endl;
//     Info << "subOuterRelTol_ = " << subOuterRelTol_ << endl;
//     Info << "subRelax_ = " << subRelax_ << endl;

    Info << endl;
    Info << "Relative tolerance = " << subRelTol_ << " / " << relTol_ << endl;
    if (checkSubRelax())
    {
        Info << "Relax factor = " << subRelax_  << " / " << relax_ << endl;
    }
    Info << endl;
}

void eddyCurrentControl::decreaseSubScale() const
{
    subScale_ /= 10.0;

    subOuterRelTol_ /= 10.0;

    Info << nl << "Outer relative tolerance = " << subOuterRelTol_ << endl << nl;

    subIter_ = 0;
}

const bool& eddyCurrentControl::subLoop() const
{
    subLoop_ = !checkSubConvergence();

    subIter_++;

    return subLoop_;
};

void eddyCurrentControl::reset()
{
    readDictDataIfModified();

    stop_ = false;
    loop_ = false;
    iter_ = 0;

    subLoop_ = false;
    subIter_ = -1;
    subRelTol_ = relTol_;
    subOuterRelTol_ = outerRelTol_;
    subRelax_ = relax_;
    subScale_ = pow(tol_,-1.0);
    subAdict_ = Adict_;
    subVdict_ = Vdict_;

    AVres_ = GREAT;
    AVrelRes_ = GREAT;
    Ares_ = GREAT;
    Vres_ = GREAT;

    oldAVres_ = 10*AVres_;
    oldAVrelRes_ = 10*AVrelRes_;
    oldAVresSet_ = false;

    if (!mesh3D_)
    {
        Vres_ = 0.0;

        subRelTol_ = 0.0;

        solDir_ = -mesh_[baseRegion_].geometricD();
    }

    // Update subDict copies
    resetSubDictionaries();

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

eddyCurrentControl::eddyCurrentControl
(
    const regionFvMesh& mesh,
    const word& name,
    const bool& master
)
:
    solverControl<regionFvMesh>
    (
        mesh,
        name,
        master,
        mesh.regionIndex(polyMesh::defaultRegion)
    ),
    conductorRegionName_
    (
        word(propDict_.subDict("regions").lookup("conductor"))
    ),
    conductorRegion_
    (
        mesh_.regionIndex(conductorRegionName_)
    ),
    interfacePatchName_
    (
        word(propDict_.subDict("patches").lookup("interface"))
    ),
    interfacePatch_
    (
        mesh_[conductorRegion_].boundaryMesh().findPatchID(interfacePatchName_)
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
    conductorSolutionDict_(mesh_[conductorRegion_].solutionDict()),
    AVdict_(baseSolutionDict_.subDict("solvers").subDict("AV")),
    Adict_(baseSolutionDict_.subDict("solvers").subDict("A")),
    Vdict_(conductorSolutionDict_.subDict("solvers").subDict("V")),
    tol_(AVdict_.lookupOrDefault<scalar>("tolerance", 1e-04)),
    relTol_(AVdict_.lookupOrDefault<scalar>("relTol", 0.5)),
    outerRelTol_(AVdict_.lookupOrDefault<scalar>("outerRelTol", 0.05)),
    outerMaxIter_(AVdict_.lookupOrDefault<int>("outerMaxIter", 100)),
    relax_(AVdict_.lookupOrDefault<scalar>("relax", 1.0)),
    mesh3D_((mesh_[baseRegion_].nGeometricD() == 3)),
    solDir_(mesh_[baseRegion_].geometricD()),
    stop_(false),
    loop_(false),
    iter_(0),
    subLoop_(false),
    subIter_(-1),
    subRelTol_(relTol_),
    subOuterRelTol_(outerRelTol_),
    subRelax_(relax_),
    subScale_(pow(tol_,-1.0)),
    subAdict_(Adict_),
    subVdict_(Vdict_),
    AVres_(GREAT),
    AVrelRes_(GREAT),
    Ares_(GREAT),
    Vres_(GREAT),
    oldAVres_(10*AVres_),
    oldAVrelRes_(10*AVrelRes_),
    oldAVresSet_(false)
{
    reset();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void eddyCurrentControl::relax(volVectorField& vf) const
{
    if (checkSubRelax())
    {
        if (!firstIteration())
        {
            Info << "Relax " << vf.name()
                << " with a factor = " << subRelax_ << endl;

            vf.relax(subRelax_);
        }

        vf.storePrevIter();
    }
}

void eddyCurrentControl::relax(const regionVolVectorField& rvf) const
{
    if (checkSubRelax())
    {
        if (!firstIteration())
        {
            Info << "Relax " << rvf.name()
                << " with a factor = " << subRelax_ << endl;

            forAll(rvf.mesh().regionNames(), regionI)
            {
                rvf[regionI].relax(subRelax_);
            }
        }

        forAll(rvf.mesh().regionNames(), regionI)
        {
            rvf[regionI].storePrevIter();
        }
    }
}

const bool& eddyCurrentControl::loop()
{
    bool converged = checkConvergence();
    bool iterBelowMax = checkMaxIterations();

//     // Update data if last subiteration has converged
//     if (!converged && subConverged)
//     {
//         readDictDataIfModified();
//     }

    // Loop if not converged and below max iterations
    loop_ = !converged && iterBelowMax;

    if (loop_)
    {
        // Reset calculation switch for old residual
        oldAVresSet_ = false;

        if (!mesh3D_)
        {
            iter_++;
        }
        else
        {
            if (!subLoop())
            {
                iter_++;

                decreaseSubScale();
            }

            decreaseSubTolerance();
        }
    }
    else
    {
        reset();
    }

    return loop_;
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

