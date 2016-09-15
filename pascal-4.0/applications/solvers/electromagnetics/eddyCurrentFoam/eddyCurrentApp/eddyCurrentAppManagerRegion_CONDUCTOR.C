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

void Foam::eddyCurrentApp::Manager::Region_CONDUCTOR::Settings::read() const
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::eddyCurrentApp::Manager::Region_CONDUCTOR::Storage::Item_ddtARePhi::create() const
{
    set
    (
        new surfaceScalarField
        (
            IOobject
            (
                "ddtARePhi",
                time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar
            (
                word(),
                dimVoltage*dimLength,
                0
            ),
            calculatedFvPatchScalarField::typeName
        )
    );
}


void Foam::eddyCurrentApp::Manager::Region_CONDUCTOR::Storage::Item_ddtAImPhi::create() const
{
    set
    (
        new surfaceScalarField
        (
            IOobject
            (
                "ddtAImPhi",
                time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar
            (
                word(),
                dimVoltage*dimLength,
                0
            ),
            calculatedFvPatchScalarField::typeName
        )
    );
}


void Foam::eddyCurrentApp::Manager::Region_CONDUCTOR::Storage::Item_VRe::create() const
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


void Foam::eddyCurrentApp::Manager::Region_CONDUCTOR::Storage::Item_VIm::create() const
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


void Foam::eddyCurrentApp::Manager::Region_CONDUCTOR::Storage::Item_jRe::create() const
{
    IOobject IOo
    (
        "jRe",
        time().timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    );

    if (!dict().lookupOrDefault<bool>("write", true))
    {
        IOo.writeOpt() = IOobject::NO_WRITE;
    }

    set
    (
        new volVectorField
        (
            IOo,
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


void Foam::eddyCurrentApp::Manager::Region_CONDUCTOR::Storage::Item_jIm::create() const
{
    IOobject IOo
    (
        "jIm",
        time().timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    );

    if (!dict().lookupOrDefault<bool>("write", true))
    {
        IOo.writeOpt() = IOobject::NO_WRITE;
    }

    set
    (
        new volVectorField
        (
            IOo,
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::eddyCurrentApp::Manager::Region_CONDUCTOR::Storage::create() const
{
    item_ddtARePhi().enable();
    item_ddtAImPhi().enable();

    item_VRe().setState(control().meshIs3D());
    item_VIm().setState(control().meshIs3D());

    item_jRe().enable();
    item_jIm().enable();
}


// ************************************************************************* //

