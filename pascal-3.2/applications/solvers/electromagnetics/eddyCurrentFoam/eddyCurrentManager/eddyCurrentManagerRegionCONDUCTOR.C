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

#include "eddyCurrentManager.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void eddyCurrentManager::Region_CONDUCTOR::Settings::read() const
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void eddyCurrentManager::Region_CONDUCTOR::Storage::Item_VRe::create
(const word& ccase) const
{
    set
    (
        new volScalarField
        (
            IOobject
            (
                "VRe",
                time().timeName(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh()
        )
    );
}


void eddyCurrentManager::Region_CONDUCTOR::Storage::Item_VIm::create
(const word& ccase) const
{
    set
    (
        new volScalarField
        (
            IOobject
            (
                "VIm",
                time().timeName(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh()
        )
    );
}


void eddyCurrentManager::Region_CONDUCTOR::Storage::Item_VReAverage::create
(const word& ccase) const
{
    set
    (
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "VReAverage",
                time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            dimensionedScalar(
                word(),
                dimVoltage,
                gAverage(storage().VRe().internalField())
            )
        )
    );
}


void eddyCurrentManager::Region_CONDUCTOR::Storage::Item_VImAverage::create
(const word& ccase) const
{
    set
    (
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "VImAverage",
                time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            dimensionedScalar(
                word(),
                dimVoltage,
                gAverage(storage().VIm().internalField())
            )
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void eddyCurrentManager::Region_CONDUCTOR::Storage::create
(const word& ccase) const
{
    item_VRe().enable();
    item_VIm().enable();
    
    item_VReAverage().enable();
    item_VImAverage().enable();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

