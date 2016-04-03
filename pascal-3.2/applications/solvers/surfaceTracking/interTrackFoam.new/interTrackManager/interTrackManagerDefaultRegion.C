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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void interTrackManager::DefaultRegion::Settings::read() const
{
// TODO: Add Settings dictionary
    const dictionary& settingsDict = this->settingsRegionDict();

    UpCoupled = settingsDict.lookupOrDefault("UpCoupled", false);
    UpDirectForce = settingsDict.lookupOrDefault("UpDirectForce", false);
}


void interTrackManager::DefaultRegion::Storage::init_g(const word& init) const
{
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


void interTrackManager::DefaultRegion::Storage::init_p(const word& init) const
{
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


void interTrackManager::DefaultRegion::Storage::init_U(const word& init) const
{
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


void interTrackManager::DefaultRegion::Storage::init_phi(const word& init) const
{
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


void interTrackManager::DefaultRegion::Storage::init_rho(const word& init) const
{
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
            dimensionedScalar(word(), dimDensity, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


void interTrackManager::DefaultRegion::Storage::init_mu(const word& init) const
{
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
            dimensionedScalar(word(), dimMass/dimLength/dimTime, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


void interTrackManager::DefaultRegion::Storage::init_F(const word& init) const
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


void interTrackManager::DefaultRegion::Storage::init_fluidIndicator(const word& init) const
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


void interTrackManager::DefaultRegion::Storage::init_transport(const word& init) const
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


void interTrackManager::DefaultRegion::Storage::init_turbulence(const word& init) const
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


void interTrackManager::DefaultRegion::Storage::init_interface(const word& init) const
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
            this->mesh(),
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


void interTrackManager::DefaultRegion::Storage::init() const
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

