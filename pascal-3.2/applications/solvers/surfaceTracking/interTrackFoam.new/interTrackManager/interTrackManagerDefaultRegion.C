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
#include "zeroGradientFvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void interTrackManager::DefaultRegion::Settings::read() const
{
    UpCoupled =
        dict().lookupOrDefault("UpCoupled", false);
    UEqnVolumeForce =
        dict().lookupOrDefault("UEqnVolumeForce", false);
    snGradpFromFlux =
        dict().lookupOrDefault("snGradpFromFlux", true);

    if
    (
        UpCoupled
    && !UEqnVolumeForce
    )
    {
        FatalErrorIn("interTrackManager::DefaultRegion::Settings::read() : ")
            << "Both coupled solution of U and p and indirect flux "
            << "representation of enabled volume force is currently "
            << "not implemented!"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void interTrackManager::DefaultRegion::Storage::init_g
(const word& init) const
{
    gPtr_.set
    (
        new uniformDimensionedVectorField
        (
            IOobject
            (
                "g",
                time().constant(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );
}


void interTrackManager::DefaultRegion::Storage::init_p
(const word& init) const
{
    pPtr_.set
    (
        new volScalarField
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
        )
    );
}


void interTrackManager::DefaultRegion::Storage::init_U
(const word& init) const
{
    UPtr_.set
    (
        new volVectorField
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
        )
    );
}


void interTrackManager::DefaultRegion::Storage::init_Up
(const word& init) const
{
    UpPtr_.set
    (
        new volVector4Field
        (
            IOobject
            (
                "Up",
                time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedVector4
            (
                word(),
                dimless,
                vector4::zero
            )
        )
    );
}


void interTrackManager::DefaultRegion::Storage::init_phi
(const word& init) const
{
    phiPtr_.set
    (
        new surfaceScalarField
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
        )
    );

    phiPtr_->oldTime();
}


void interTrackManager::DefaultRegion::Storage::init_rho
(const word& init) const
{
    rhoPtr_.set
    (
        new volScalarField
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
            dimensionedScalar
            (
                word(),
                dimDensity,
                0
            ),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


void interTrackManager::DefaultRegion::Storage::init_mu
(const word& init) const
{
    muPtr_.set
    (
        new volScalarField
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
            dimensionedScalar
            (
                word(),
                dimMass/dimLength/dimTime,
                0
            ),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


void interTrackManager::DefaultRegion::Storage::init_F
(const word& init) const
{
    if (init == "default")
    {
        FPtr_.set
        (
            new volVectorField
            (
                IOobject
                (
                    "F",
                    time().timeName(),
                    mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh()
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
                    time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimensionedVector
                (
                    "F",
                    dimForce/dimVolume,
                    vector::zero
                ),
                calculatedFvPatchVectorField::typeName
            )
        );
    }
}


void interTrackManager::DefaultRegion::Storage::init_fluidIndicator
(const word& init) const
{
    fluidIndicatorPtr_.set
    (
        new volScalarField
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
        )
    );
}


void interTrackManager::DefaultRegion::Storage::init_transport
(const word& init) const
{
    transportPtr_.set
    (
        new twoPhaseMixture
        (
            U(),
            phi(),
            fluidIndicator().name()
        )
    );
}


void interTrackManager::DefaultRegion::Storage::init_turbulence
(const word& init) const
{
    turbulencePtr_ =
    (
        incompressible::turbulenceModel::New
        (
            U(),
            phi(),
            transport()
        )
    );
}


void interTrackManager::DefaultRegion::Storage::init_interface
(const word& init) const
{
    word interfacePrefix;

    if (!args().optionReadIfPresent("prefix", interfacePrefix))
    {
        interfacePrefix = trackedSurface::typeName;
    }

    // Include turbulence only if it is activated
    incompressible::turbulenceModel* turbPtr = NULL;
    if(is_turbulence())
    {
        turbPtr = turbulencePtr();
    }

// TODO: Add more constructors and simplify
//       Make it use real rho and mu fields!
    interfacePtr_.set
    (
        new trackedSurface
        (
            mesh(),
            rho(),
            U(),
            p(),
            phi(),
            NULL,
            gPtr(),
            transportPtr(),
            turbPtr,
            NULL,
            interfacePrefix
        )
    );
}

// TODO: makeOptional_ has to be connected to settings!
void interTrackManager::DefaultRegion::Storage::create() const
{
    make_g();
    make_p();
    make_U();
    make_phi();
    make_rho();
    make_mu();
    makeOptional_F();
    make_fluidIndicator();
    make_transport();

    // Init density from transport (viscosity) model
    {
        dimensionedScalar rho1 = transport().rho1();
        dimensionedScalar rho2 = transport().rho2();

        rho() = fluidIndicator()*(rho1 - rho2) + rho2;
        rho().correctBoundaryConditions();
    }

    // Init (dynamic) viscosity from transport (viscosity) model and density
    {
        dimensionedScalar nu1 = dimensionedScalar
        (
            transport().nuModel1().viscosityProperties().lookup("nu")
        );
        dimensionedScalar nu2 = dimensionedScalar
        (
            transport().nuModel2().viscosityProperties().lookup("nu")
        );

        mu() = (fluidIndicator()*(nu1 -nu2) + nu2) * rho();
        mu().correctBoundaryConditions();
    }

    makeOptional_turbulence();
    make_interface();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

