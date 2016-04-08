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

void eddyCurrentManager::DefaultRegion::Settings::read() const
{
// TODO
//     manager().settings();
//     region().settings();
//
//     manager().storage().f0();
//     manager().regions().defaultRegion().storage().j0Re();
//     region().storage().j0Re();
//     manager().regions().conductorRegion().storage().VRe();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void eddyCurrentManager::DefaultRegion::Storage::Item_j0Re::create
(const word& init) const
{
    set
    (
        new volVectorField
        (
            IOobject
            (
                "j0Re",
                time().timeName(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedVector
            (
                word(),
                dimCurrent/dimArea,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName
        )
    );
}


void eddyCurrentManager::DefaultRegion::Storage::Item_j0Im::create
(const word& init) const
{
    set
    (
        new volVectorField
        (
            IOobject
            (
                "j0Im",
                time().timeName(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh()
        )
    );
}


void eddyCurrentManager::DefaultRegion::Storage::Item_jRe::create
(const word& init) const
{
    set
    (
        new volVectorField
        (
            IOobject
            (
                "jRe",
                time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedVector
            (
                word(),
                dimCurrent/dimArea,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName
        )
    );
}


void eddyCurrentManager::DefaultRegion::Storage::Item_jIm::create
(const word& init) const
{
    set
    (
        new volVectorField
        (
            IOobject
            (
                "jIm",
                time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedVector
            (
                word(),
                dimCurrent/dimArea,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName
        )
    );
}


void eddyCurrentManager::DefaultRegion::Storage::Item_BRe::create
(const word& init) const
{
    set
    (
        new volVectorField
        (
            IOobject
            (
                "BRe",
                time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedVector
            (
                word(),
                dimVoltage*dimTime/dimArea,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName
        )
    );
}


void eddyCurrentManager::DefaultRegion::Storage::Item_BIm::create
(const word& init) const
{
    set
    (
        new volVectorField
        (
            IOobject
            (
                "BIm",
                time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedVector
            (
                word(),
                dimVoltage*dimTime/dimArea,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName
        )
    );
}


void eddyCurrentManager::DefaultRegion::Storage::create() const
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

