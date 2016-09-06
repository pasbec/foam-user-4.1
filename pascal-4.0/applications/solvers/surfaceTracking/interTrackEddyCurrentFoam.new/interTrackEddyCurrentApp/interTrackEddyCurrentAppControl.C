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

#include "interTrackEddyCurrentAppControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::interTrackEddyCurrentApp::Control, 0);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::interTrackEddyCurrentApp::Control::read()
{
    // Read settings
    emSettings_.enabled =
        emUpdateSettingsDict_.lookupOrDefault<bool>("enabled", true);
    emSettings_.outputTimeIndexCycle =
        emUpdateSettingsDict_.lookupOrDefault<int>("outputTimeIndexCycle", 2147483647);
    emSettings_.timeIndexCycle =
        emUpdateSettingsDict_.lookupOrDefault<int>("timeIndexCycle", 2147483647);
    emSettings_.timeCycle =
        emUpdateSettingsDict_.lookupOrDefault<scalar>("timeCycle", VGREAT);
    emSettings_.relDeltaAmax =
        emUpdateSettingsDict_.lookupOrDefault<scalar>("relDeltaAmax", 0.01);

    if (debug > 1)
    {
        Info << "interTrackEddyCurrentApp::Control::read() : "
            << "emUpdateSettingsDict = " << emUpdateSettingsDict_;
    }

    // Read data
    emUpdateDataDict_.readIfModified();
    emUpdateData_.update = false;
    emUpdateData_.counter =
        emUpdateDataDict_.lookupOrAddDefault<int>("counter", 0);
    emUpdateData_.outputTimeIndex =
        emUpdateDataDict_.lookupOrAddDefault<int>("outputTimeIndex", 0);
    emUpdateData_.lastTime =
        emUpdateDataDict_.lookupOrAddDefault<scalar>("lastTime", -VGREAT);

    if (debug > 1)
    {
        Info << "interTrackEddyCurrentApp::Control::read() : "
            << "emUpdateDataDict = " << emUpdateDataDict_;
    }
}

bool Foam::interTrackEddyCurrentApp::Control::criteriaSatisfied()
{
    emUpdateData_.update = false;

    if (emSettings_.enabled)
    {
        // Increase output index if due and write to dictionary
        if (mesh_.time().outputTime())
        {
            emUpdateData_.outputTimeIndex += 1;
            emUpdateDataDict_.set<int>
            (
                "outputTimeIndex",
                emUpdateData_.outputTimeIndex
            );
        }

        emUpdateData_.update = update();

        if (emUpdateData_.update)
        {
            // Increase update counter and write to dictionary
            emUpdateData_.counter += 1;
            emUpdateDataDict_.set<int>
            (
                "counter",
                emUpdateData_.counter
            );

            // Update last update time and write to dictionary
            emUpdateData_.lastTime = mesh_.time().value();
            emUpdateDataDict_.set<scalar>
            (
                "lastTime",
                emUpdateData_.lastTime
            );
        }
    }

    return !emUpdateData_.update;
}

bool Foam::interTrackEddyCurrentApp::Control::updateRelDeltaA() const
{
    regionVolVectorField C
    (
        IOobject
        (
            "C",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector
        (
            word(),
            dimLength,
            vector::zero
        ),
        calculatedFvPatchVectorField::typeName
    );

    C[Region::DEFAULT] = mesh_[Region::DEFAULT].C();
    C[Region::FLUID] = mesh_[Region::FLUID].C();
    C.rmapInteralField(Region::FLUID);
    C.mapCopyInternal(Region::CONDUCTOR);

    scalarField magSqrDeltaC =
        magSqr
        (
            C[Region::CONDUCTOR].internalField()
          - prevC_.internalField()
        );

    const volVectorField& jRe =
        mesh_[Region::CONDUCTOR].lookupObject<volVectorField> ("jRe");
    const volVectorField& jIm =
        mesh_[Region::CONDUCTOR].lookupObject<volVectorField> ("jIm");

    scalarField magSqrj = magSqr(jRe) + magSqr(jIm);

    const volVectorField& ARe =
        mesh_[Region::CONDUCTOR].lookupObject<volVectorField> ("ARe");
    const volVectorField& AIm =
        mesh_[Region::CONDUCTOR].lookupObject<volVectorField> ("AIm");

    scalarField magSqrA(mesh_[Region::CONDUCTOR].C().size(),0.0);

    if
    (
        mesh_[Region::CONDUCTOR].foundObject<volVectorField>("A0Re")
     && mesh_[Region::CONDUCTOR].foundObject<volVectorField>("A0Im")
    )
    {

        const volVectorField& A0Re =
            mesh_[Region::CONDUCTOR].lookupObject<volVectorField> ("A0Re");
        const volVectorField& A0Im =
            mesh_[Region::CONDUCTOR].lookupObject<volVectorField> ("A0Im");

        magSqrA = magSqr(ARe+A0Re) + magSqr(AIm+A0Im);
    }
    else
    {
        magSqrA = magSqr(ARe) + magSqr(AIm);
    }

    relDeltaA_.internalField() =
        physicalConstant::mu0.value()
        * magSqrDeltaC * sqrt(magSqrj / (magSqrA + SMALL));
    relDeltaA_.correctBoundaryConditions();

    if (gMax(relDeltaA_) > emSettings_.relDeltaAmax)
    {
        prevC_ = mesh_[Region::CONDUCTOR].C();
        prevC_.correctBoundaryConditions();

        return true;
    }

    return false;
}

bool Foam::interTrackEddyCurrentApp::Control::update() const
{
    if (debug)
    {
        Info << "interTrackEddyCurrentApp::Control::update() : "
            << "updateZeroCounter() = " << updateZeroCounter() << endl;
        Info << "interTrackEddyCurrentApp::Control::update() : "
            << "updateOutputTimeIndex() = " << updateOutputTimeIndex() << endl;
        Info << "interTrackEddyCurrentApp::Control::update() : "
            << "updateTimeIndex() = " << updateTimeIndex() << endl;
        Info << "interTrackEddyCurrentApp::Control::update() : "
            << "updateTime() = " << updateTime() << endl;
        Info << "interTrackEddyCurrentApp::Control::update() : "
            << "updateRelDeltaA() = " << updateRelDeltaA() << endl;
    }

    return updateZeroCounter()
        || updateOutputTimeIndex()
        || updateTimeIndex()
        || updateTime()
        || updateRelDeltaA();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interTrackEddyCurrentApp::Control::Control
(
    regionDynamicFvMesh& mesh,
    const dictionary& propertiesDict,
    const word& dictName
)
:
    solutionControl(mesh[Region::DEFAULT], dictName),
    mesh_(mesh),
    prevC_
    (
        const_cast<volVectorField&>
        (
            mesh_[Region::CONDUCTOR].lookupObject<volVectorField> ("emPrevC")
        )
    ),
    relDeltaA_
    (
        const_cast<volScalarField&>
        (
            mesh_[Region::CONDUCTOR].lookupObject<volScalarField> ("emRelDeltaA")
        )
    ),
    propertiesDict_
    (
        propertiesDict
    ),
    emUpdateSettingsDict_
    (
        propertiesDict_.subDict("emUpdate")
    ),
    emUpdateDataDict_
    (
        IOobject
        (
            "emUpdate",
            mesh.time().timeName(),
            "uniform",
            mesh[Region::DEFAULT],
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE,
            true
        )
    )
{
    read();

    prevC_ = mesh_[Region::CONDUCTOR].C();
    prevC_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interTrackEddyCurrentApp::Control::~Control()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::interTrackEddyCurrentApp::Control::loop()
{
    read();

    return !criteriaSatisfied();
}

// ************************************************************************* //
