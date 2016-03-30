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

#include "interTrackManager.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void interTrackManager::storage::create_g(const word init) const
{
    solverManagerStorage_info(g);

    gPtr_.set
    (
        new uniformDimensionedVectorField
        (
            IOobject
            (
                "g",
                this->time().constant(),
                this->mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );
}


void interTrackManager::storage::create_p(const word init) const
{
    solverManagerStorage_info(p);

    pPtr_.set
    (
        new volScalarField
        (
            IOobject
            (
                "p",
                this->time().timeName(),
                this->mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            this->mesh()
        )
    );
}


void interTrackManager::storage::create_U(const word init) const
{
    solverManagerStorage_info(U);

    UPtr_.set
    (
        new volVectorField
        (
            IOobject
            (
                "U",
                this->time().timeName(),
                this->mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            this->mesh()
        )
    );
}


void interTrackManager::storage::create_phi(const word init) const
{
    solverManagerStorage_info(phi);

    solverManagerStorage_assert(phi, U);

    phiPtr_.set
    (
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                this->time().timeName(),
                this->mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(U()) & this->mesh().Sf()
        )
    );

    phiPtr_->oldTime();
}


void interTrackManager::storage::create_rho(const word init) const
{
    solverManagerStorage_info(rho);

    rhoPtr_.set
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(word(), dimDensity, 0)
        )
    );
}


void interTrackManager::storage::create_mu(const word init) const
{
    solverManagerStorage_info(mu);

    muPtr_.set
    (
        new volScalarField
        (
            IOobject
            (
                "mu",
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(word(), dimMass/dimLength/dimTime, 0)
        )
    );
}


void interTrackManager::storage::create_F(const word init) const
{
    solverManagerStorage_info(F);

    if (init == "default")
    {
        FPtr_.set
        (
            new volVectorField
            (
                IOobject
                (
                    "F",
                    this->time().timeName(),
                    this->mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                this->mesh()
            )
        );
    }
    else if (init == "calculated")
    {
        FPtr_.set
        (
            new volVectorField
            (
                IOobject
                (
                    "F",
                    this->time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedVector
                (
                    "F",
                    dimMass/pow(dimLength,2)/pow(dimTime,2),
                    pTraits<vector>::zero
                ),
                calculatedFvPatchVectorField::typeName
            )
        );
    }
}


void interTrackManager::storage::create_fluidIndicator(const word init) const
{
    solverManagerStorage_info(fluidIndicator);

    solverManagerStorage_assert(fluidIndicator, interface);

    if(args().options().found("parallel"))
    {
        fluidIndicatorPtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "fluidIndicator",
                    this->time().timeName(),
                    this->mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                this->mesh()
            )
        );
    }
    else
    {
        fluidIndicatorPtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "fluidIndicator",
                    this->time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                interface().fluidIndicator()
            )
        );
    }

    fluidIndicatorPtr_->correctBoundaryConditions();
    fluidIndicatorPtr_->write();
}


void interTrackManager::storage::create_interface(const word init) const
{
    solverManagerStorage_info(interface);

    solverManagerStorage_assert(interface, p);
    solverManagerStorage_assert(interface, U);
    solverManagerStorage_assert(interface, phi);
    solverManagerStorage_assert(interface, rho);

    word interfacePrefix;

    if (!args().optionReadIfPresent("prefix", interfacePrefix))
    {
        interfacePrefix = trackedSurface::typeName;
    }

// TODO: Add more constructors and simplify
    interfacePtr_.set
    (
        new trackedSurface
        (
            this->mesh(),
            rho(),
            U(),
            p(),
            phi(),
            NULL,
            gPtr(),
            NULL,
            NULL,
            NULL,
            interfacePrefix
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void interTrackManager::storage::init() const
{
    g();
    p();
    U();
    phi();
    rho();
    mu();
    F("calculated");

    // To create the interface we need rho first
    interface();

    // The fluid indicator may be created from interface
    fluidIndicator();

    // The density rho is derived from the fluid indicator
    // and the interface density values
    {
        rho() = fluidIndicator()
           *(
                interface().rhoFluidA()
              - interface().rhoFluidB()
            )
          + interface().rhoFluidB();

        rho().correctBoundaryConditions();
    }

    // The viscosity mu is derived from the fluid indicator
    // and the interface viscosity values
    {
        mu() = fluidIndicator()
           *(
                interface().muFluidA()
              - interface().muFluidB()
            )
          + interface().muFluidB();

        mu().correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

