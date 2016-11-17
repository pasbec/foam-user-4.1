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

#include "buoyantBoussinesqPimpleEddyCurrentAppManager.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager::Region_FLUID::Settings::read() const
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager::Region_FLUID::Storage::Item_g::create() const
{
    set
    (
        new uniformDimensionedVectorField
        (
            IOobject
            (
                name(),
                time().constant(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );
}


void Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager::Region_FLUID::Storage::Item_beta::create() const
{
    set
    (
        new uniformDimensionedScalarField
        (
            IOobject
            (
                name(),
                time().constant(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            dimensionedScalar
            (
                manager().pimpleAppManager().regions().region_DEFAULT().storage().transport().lookup("beta")
            )
        )
    );
}


void Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager::Region_FLUID::Storage::Item_TRef::create() const
{
    set
    (
        new uniformDimensionedScalarField
        (
            IOobject
            (
                name(),
                time().constant(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            dimensionedScalar
            (
                manager().pimpleAppManager().regions().region_DEFAULT().storage().transport().lookup("TRef")
            )
        )
    );
}


void Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager::Region_FLUID::Storage::Item_Pr::create() const
{
    set
    (
        new uniformDimensionedScalarField
        (
            IOobject
            (
                name(),
                time().constant(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            dimensionedScalar
            (
                manager().pimpleAppManager().regions().region_DEFAULT().storage().transport().lookup("Pr")
            )
        )
    );
}


void Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager::Region_FLUID::Storage::Item_Prt::create() const
{
    set
    (
        new uniformDimensionedScalarField
        (
            IOobject
            (
                name(),
                time().constant(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            dimensionedScalar
            (
                manager().pimpleAppManager().regions().region_DEFAULT().storage().transport().lookup("Prt")
            )
        )
    );
}


void Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager::Region_FLUID::Storage::Item_kappat::create() const
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


// void Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager::Region_FLUID::Storage::Item_rhok::create() const
// {
//     Region_THERMAL& region_THERMAL = manager().regions().region_THERMAL();
//
//     set
//     (
//         new volScalarField
//         (
//             IOobject
//             (
//                 name(),
//                 time().timeName(),
//                 mesh(),
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE
//             ),
//             1.0 - storage().beta()*(Region_THERMAL().storage().T() - storage().TRef()),
//             calculatedFvPatchScalarField::typeName
//         )
//     );
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager::Region_FLUID::Storage::create() const
{
    item_g().enable();
    item_beta().enable();
    item_TRef().enable();
    item_Pr().enable();
    item_Prt().enable();
    item_kappat().enable();
//     item_rhok().enable();
}


// ************************************************************************* //

