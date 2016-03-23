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

#include "interTrackControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void interTrackControl::storage::create_p()
{
    solverControlStorage_info(p);

    pPtr_ = new volScalarField
    (
        IOobject
        (
            "p",
            time().timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    );
}


void interTrackControl::storage::create_U()
{
    solverControlStorage_info(U);

    UPtr_ = new volVectorField
    (
        IOobject
        (
            "U",
            time().timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    );
}


void interTrackControl::storage::create_phi()
{
    solverControlStorage_info(phi);

    solverControlStorage_assert(phi, U);

    phiPtr_ = new surfaceScalarField
    (
        IOobject
        (
            "phi",
            time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(U()) & mesh().Sf()
    );

    phiPtr_->oldTime();
}


void interTrackControl::storage::create_rho()
{
    solverControlStorage_info(rho);

    rhoPtr_ = new volScalarField
    (
        IOobject
        (
            "rho",
            time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar(word(), dimDensity, 0)
    );
}


void interTrackControl::storage::create_mu()
{
    solverControlStorage_info(mu);

    muPtr_ = new volScalarField
    (
        IOobject
        (
            "mu",
            time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar(word(), dimMass/dimLength/dimTime, 0)
    );
}


void interTrackControl::storage::create_fluidIndicator()
{
    solverControlStorage_info(fluidIndicator);

    solverControlStorage_assert(fluidIndicator, interface);

    if(args().options().found("parallel"))
    {
        fluidIndicatorPtr_ = new volScalarField
        (
            IOobject
            (
                "fluidIndicator",
                time().timeName(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh()
        );
    }
    else
    {
        fluidIndicatorPtr_ = new volScalarField
        (
            IOobject
            (
                "fluidIndicator",
                time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            interface().fluidIndicator()
        );
    }

    fluidIndicatorPtr_->correctBoundaryConditions();
    fluidIndicatorPtr_->write();
}


void interTrackControl::storage::create_interface()
{
    solverControlStorage_info(interface);

    solverControlStorage_assert(interface, p);
    solverControlStorage_assert(interface, U);
    solverControlStorage_assert(interface, phi);
    solverControlStorage_assert(interface, rho);

    word interfacePrefix;

    if (!args().optionReadIfPresent("prefix", interfacePrefix))
    {
        interfacePrefix = trackedSurface::typeName;
    }

    interfacePtr_ = new trackedSurface
    (
        mesh(),
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

void interTrackControl::storage::init()
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

interTrackControl::storage::storage(interTrackControl& control)
:
    solverControl<dynamicFvMesh>::storage
    (
        control.args(), control.time(), control.mesh()
    ),
    control_(control),
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

