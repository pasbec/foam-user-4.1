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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(eddyCurrentManager, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void eddyCurrentManager::Settings::read() const
{
    debug =
        dict().lookupOrDefault
        (
            "debug", Switch(eddyCurrentManager::debug)
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void eddyCurrentManager::Storage::init_f0
(const word& init) const
{
    f0Ptr_.set
    (
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "f0",
                time().constant(),
                time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );
}


void eddyCurrentManager::Storage::init_omega0
(const word& init) const
{
    omega0Ptr_.set
    (
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "omega0",
                time().constant(),
                time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mathematicalConstant::twoPi*f0()
        )
    );
}


void eddyCurrentManager::Storage::init_sigma
(const word& init) const
{
    sigmaPtr_.set
    (
        new regionVolScalarField
        (
            IOobject
            (
                "sigma",
                time().timeName(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh()
        )
    );
}


void eddyCurrentManager::Storage::init_ARe
(const word& init) const
{
    ARePtr_.set
    (
        new regionVolVectorField
        (
            IOobject
            (
                "ARe",
                time().timeName(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh()
        )
    );
}


void eddyCurrentManager::Storage::init_AIm
(const word& init) const
{
    AImPtr_.set
    (
        new regionVolVectorField
        (
            IOobject
            (
                "AIm",
                time().timeName(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh()
        )
    );
}


void eddyCurrentManager::Storage::init_VReGrad
(const word& init) const
{
    VReGradPtr_.set
    (
        new regionVolVectorField
        (
            IOobject
            (
                "VReGrad",
                time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedVector
            (
                word(),
                dimVoltage/dimLength,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName
        )
    );
}


void eddyCurrentManager::Storage::init_VImGrad
(const word& init) const
{
    VImGradPtr_.set
    (
        new regionVolVectorField
        (
            IOobject
            (
                "VImGrad",
                time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedVector
            (
                word(),
                dimVoltage/dimLength,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName
        )
    );
}


void eddyCurrentManager::Storage::init_FL
(const word& init) const
{
    FLPtr_.set
    (
        new regionVolVectorField
        (
            IOobject
            (
                "FL",
                time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedVector
            (
                word(),
                dimForce/dimVolume,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName
        )
    );
}


void eddyCurrentManager::Storage::init_pB
(const word& init) const
{
    pBPtr_.set
    (
        new regionVolScalarField
        (
            IOobject
            (
                "pB",
                time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar
            (
                word(),
                dimPressure,
                0
            ),
            calculatedFvPatchScalarField::typeName
        )
    );
}


void eddyCurrentManager::Storage::init() const
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void eddyCurrentManager::next() const
{}


void eddyCurrentManager::write() const
{}


void eddyCurrentManager::finalize() const
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

eddyCurrentManager::eddyCurrentManager
(
    const argList& args,
    Time& time,
    regionFvMesh& mesh,
    const word& name,
    const bool& master
)
:
    solverManagerCore<regionFvMesh>
    (
        args, time, mesh, name, master
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

