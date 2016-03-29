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

void interTrackManager::storage::create_p()
{
    solverManagerStorage_info(p);

    pPtr_ = new volScalarField
    (
        IOobject
        (
            "p",
            manager().time().timeName(),
            manager().mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        manager().mesh()
    );
}


void interTrackManager::storage::create_U()
{
    solverManagerStorage_info(U);

    UPtr_ = new volVectorField
    (
        IOobject
        (
            "U",
            manager().time().timeName(),
            manager().mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        manager().mesh()
    );
}


void interTrackManager::storage::create_phi()
{
    solverManagerStorage_info(phi);

    solverManagerStorage_assert(phi, U);

    phiPtr_ = new surfaceScalarField
    (
        IOobject
        (
            "phi",
            manager().time().timeName(),
            manager().mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(U()) & manager().mesh().Sf()
    );

    phiPtr_->oldTime();
}


void interTrackManager::storage::create_rho()
{
    solverManagerStorage_info(rho);

    rhoPtr_ = new volScalarField
    (
        IOobject
        (
            "rho",
            manager().time().timeName(),
            manager().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        manager().mesh(),
        dimensionedScalar(word(), dimDensity, 0)
    );
}


void interTrackManager::storage::create_mu()
{
    solverManagerStorage_info(mu);

    muPtr_ = new volScalarField
    (
        IOobject
        (
            "mu",
            manager().time().timeName(),
            manager().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        manager().mesh(),
        dimensionedScalar(word(), dimMass/dimLength/dimTime, 0)
    );
}


void interTrackManager::storage::create_fluidIndicator()
{
    solverManagerStorage_info(fluidIndicator);

    solverManagerStorage_assert(fluidIndicator, interface);

    if(args().options().found("parallel"))
    {
        fluidIndicatorPtr_ = new volScalarField
        (
            IOobject
            (
                "fluidIndicator",
                manager().time().timeName(),
                manager().mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            manager().mesh()
        );
    }
    else
    {
        fluidIndicatorPtr_ = new volScalarField
        (
            IOobject
            (
                "fluidIndicator",
                manager().time().timeName(),
                manager().mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            interface().fluidIndicator()
        );
    }

    fluidIndicatorPtr_->correctBoundaryConditions();
    fluidIndicatorPtr_->write();
}


void interTrackManager::storage::create_interface()
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

    interfacePtr_ = new trackedSurface
    (
        manager().mesh(),
        rho(),
        U(),
        p(),
        phi(),
        NULL,
        NULL,
        NULL,
        NULL,
        NULL,
        interfacePrefix
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void interTrackManager::storage::init()
{
    p();
    U();
    phi();
    rho();

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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

interTrackManager::storage::storage(interTrackManager& manager)
:
    solverManager<dynamicFvMesh>::storage
    (
        manager.args(), manager.time(), manager.mesh()
    ),
    manager_(manager),
    pPtr_(NULL),
    UPtr_(NULL),
    phiPtr_(NULL),
    rhoPtr_(NULL),
    muPtr_(NULL),
    fluidIndicatorPtr_(NULL),
    interfacePtr_(NULL)
{
    init();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

