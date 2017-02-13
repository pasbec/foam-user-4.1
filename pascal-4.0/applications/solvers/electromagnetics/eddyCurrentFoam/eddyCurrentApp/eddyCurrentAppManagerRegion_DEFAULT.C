/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

#include "eddyCurrentAppManager.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::eddyCurrentApp::Manager::Region_DEFAULT::Settings::read() const
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::eddyCurrentApp::Manager::Region_DEFAULT::Storage::Item_gradAReNPhi::create() const
{
    set
    (
        new surfaceVectorField
        (
            IOobject
            (
                name(),
                time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedVector
            (
                word(),
                dimVoltage*dimTime,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName
        )
    );
}


void Foam::eddyCurrentApp::Manager::Region_DEFAULT::Storage::Item_gradAImNPhi::create() const
{
    set
    (
        new surfaceVectorField
        (
            IOobject
            (
                name(),
                time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedVector
            (
                word(),
                dimVoltage*dimTime,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName
        )
    );
}


void Foam::eddyCurrentApp::Manager::Region_DEFAULT::Storage::Item_f0::create() const
{
    set
    (
        new uniformDimensionedScalarField
        (
            IOobject
            (
                name(),
                time().constant(),
                time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );
}


void Foam::eddyCurrentApp::Manager::Region_DEFAULT::Storage::Item_omega0::create() const
{
    set
    (
        new uniformDimensionedScalarField
        (
            IOobject
            (
                name(),
                time().constant(),
                time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mathematicalConstant::twoPi*storage().f0()
        )
    );
}


void Foam::eddyCurrentApp::Manager::Region_DEFAULT::Storage::Item_j0Re::create() const
{
    set
    (
        new volVectorField
        (
            IOobject
            (
                name(),
                time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
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


void Foam::eddyCurrentApp::Manager::Region_DEFAULT::Storage::Item_j0Im::create() const
{
    set
    (
        new volVectorField
        (
            IOobject
            (
                name(),
                time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
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


void Foam::eddyCurrentApp::Manager::Region_DEFAULT::Storage::Item_GRe::create() const
{
    set
    (
        new volScalarField
        (
            IOobject
            (
                name(),
                time().timeName(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh()
        )
    );
}


void Foam::eddyCurrentApp::Manager::Region_DEFAULT::Storage::Item_GIm::create() const
{
    set
    (
        new volScalarField
        (
            IOobject
            (
                name(),
                time().timeName(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh()
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::eddyCurrentApp::Manager::Region_DEFAULT::Storage::create() const
{
    item_gradAReNPhi().enable();
    item_gradAImNPhi().enable();

    item_f0().enable();
    item_omega0().enable();

    item_j0Re().setState(!globalSettings().biotSavart);
    item_j0Im().setState(!globalSettings().biotSavart);

    item_GRe().enable();
    item_GIm().enable();
}


// ************************************************************************* //

