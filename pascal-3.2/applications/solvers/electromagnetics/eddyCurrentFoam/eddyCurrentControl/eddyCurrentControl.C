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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

eddyCurrentControl::eddyCurrentControl(const regionFvControl& rfvc)
:
    regIOobject
    (
        IOobject
        (
            "eddyCurrentControl",
            rfvc.rmesh().time().constant(),
            rfvc.rmesh().time().db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    control_(rfvc),
    rmesh_(rfvc.rmesh()),
    propDict_
    (
        IOdictionary
        (
            IOobject
            (
                "eddyCurrentProperties",
                rfvc.rmesh().time().constant(),
                rfvc.rmesh().time().db(),
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
    AVdict_(rmesh_[baseRegion_].solutionDict().subDict("solvers").subDict("AV")),
    Adict_(rmesh_[baseRegion_].solutionDict().subDict("solvers").subDict("A")),
    Vdict_(rmesh_[conductorRegion_].solutionDict().subDict("solvers").subDict("V")),
    minIter_(AVdict_.lookupOrDefault<int>("minIter", 0)),
    maxIter_(AVdict_.lookupOrDefault<int>("maxIter", 100)),
    tol_(AVdict_.lookupOrDefault<scalar>("tol", 1e-04))
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool eddyCurrentControl::writeData(Ostream&) const
{
    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

