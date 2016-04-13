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

#include "interTrackAppManager.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::interTrackApp::Manager::Region_DEFAULT::Settings::read() const
{
    volumeForce =
        dict().lookupOrDefault("volumeForce", false);
    turbulence =
        dict().lookupOrDefault("turbulence", false);

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
        FatalErrorIn
        (
            "Foam::interTrackApp::Manager::"
          + "Region_DEFAULT::Settings::read() : "
        )
            << "Both coupled solution of U and p and indirect flux "
            << "representation of enabled volume force is currently "
            << "not implemented!"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_g::create
(const word& ccase) const
{
    set
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


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_p::create
(const word& ccase) const
{
    set
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


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_U::create
(const word& ccase) const
{
    set
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


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_Up::create
(const word& ccase) const
{
    set
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


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_phi::create
(const word& ccase) const
{
    set
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
            linearInterpolate(storage().U()) & mesh().Sf()
        )
    );

    get()->oldTime();
}


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_rho::create
(const word& ccase) const
{
    set
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


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_mu::create
(const word& ccase) const
{
    set
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


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_F::create
(const word& ccase) const
{
    if (ccase == "default")
    {
        set
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
    else if (ccase == "calculated")
    {
        set
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


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_fluidIndicator::create
(const word& ccase) const
{
    set
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


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_transport::create
(const word& ccase) const
{
    set
    (
        new twoPhaseMixture
        (
            storage().U(),
            storage().phi(),
            storage().fluidIndicator().name()
        )
    );
}


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_turbulence::create
(const word& ccase) const
{
    set
    (
        incompressible::turbulenceModel::New
        (
            storage().U(),
            storage().phi(),
            storage().transport()
        )
    );
}


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_interface::create
(const word& ccase) const
{
    word interfacePrefix;

    if (!args().optionReadIfPresent("prefix", interfacePrefix))
    {
        interfacePrefix = trackedSurface::typeName;
    }

// TODO: Add more constructors and simplify
//       Make it use real rho and mu fields!
    set
    (
        new trackedSurface
        (
            storage().mesh(),
            storage().rho(),
            storage().U(),
            storage().p(),
            storage().phi(),
            NULL,
            storage().item_g().rawPtr(),
            storage().item_transport().rawPtr(),
            storage().item_turbulence().rawPtr(),
            NULL,
            interfacePrefix
        )
    );
}


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::create
(const word& ccase) const
{
    item_g().enable();
    item_p().enable();
    item_U().enable();
    item_phi().enable();
    item_rho().enable();
    item_mu().enable();
    item_F().setState(settings().volumeForce);

    item_fluidIndicator().enable();
    item_transport().enable();

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

    item_turbulence().setState(settings().turbulence);
    item_interface().enable();
}


// ************************************************************************* //

