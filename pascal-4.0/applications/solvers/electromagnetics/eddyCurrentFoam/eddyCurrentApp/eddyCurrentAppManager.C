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

#include "eddyCurrentAppManager.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::eddyCurrentApp::Manager, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::eddyCurrentApp::Manager::Settings::read() const
{
    Foam::eddyCurrentApp::Manager::debug =
        dict().lookupOrDefault
        (
            "debug",
            Foam::eddyCurrentApp::Manager::debug()
        );

    Foam::eddyCurrentApp::Control::debug =
        dict().lookupOrDefault
        (
            "debug",
            Foam::eddyCurrentApp::Control::debug()
        );

    biotSavart = dict().lookupOrDefault("biotSavart", false);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::eddyCurrentApp::Manager::Storage::Item_f0::create
(const word& init) const
{
    set
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


void Foam::eddyCurrentApp::Manager::Storage::Item_omega0::create
(const word& init) const
{
    set
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
            mathematicalConstant::twoPi*storage().f0()
        )
    );
}


void Foam::eddyCurrentApp::Manager::Storage::Item_sigma::create
(const word& init) const
{
    set
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


void Foam::eddyCurrentApp::Manager::Storage::Item_A0Re::create
(const word& ccase) const
{
    set
    (
        new regionVolVectorField
        (
            IOobject
            (
                "A0Re",
                time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedVector
            (
                word(),
                dimVoltage*dimTime/dimLength,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName
        )
    );
}


void Foam::eddyCurrentApp::Manager::Storage::Item_A0Im::create
(const word& ccase) const
{
    set
    (
        new regionVolVectorField
        (
            IOobject
            (
                "A0Im",
                time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedVector
            (
                word(),
                dimVoltage*dimTime/dimLength,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName
        )
    );
}


void Foam::eddyCurrentApp::Manager::Storage::Item_ARe::create
(const word& init) const
{
    set
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


void Foam::eddyCurrentApp::Manager::Storage::Item_AIm::create
(const word& init) const
{
    set
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


void Foam::eddyCurrentApp::Manager::Storage::Item_VReGrad::create
(const word& init) const
{
    set
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


void Foam::eddyCurrentApp::Manager::Storage::Item_VImGrad::create
(const word& init) const
{
    set
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


void Foam::eddyCurrentApp::Manager::Storage::Item_BRe::create
(const word& ccase) const
{
    set
    (
        new regionVolVectorField
        (
            IOobject
            (
                "BRe",
                time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedVector
            (
                word(),
                dimVoltage*dimTime/dimArea,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName
        )
    );
}


void Foam::eddyCurrentApp::Manager::Storage::Item_BIm::create
(const word& ccase) const
{
    set
    (
        new regionVolVectorField
        (
            IOobject
            (
                "BIm",
                time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedVector
            (
                word(),
                dimVoltage*dimTime/dimArea,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName
        )
    );
}


void Foam::eddyCurrentApp::Manager::Storage::Item_FL::create
(const word& init) const
{
    set
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


void Foam::eddyCurrentApp::Manager::Storage::Item_pB::create
(const word& init) const
{
    set
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


void Foam::eddyCurrentApp::Manager::Storage::create(const word& ccase) const
{
    item_f0().enable();
    item_omega0().enable();

    item_sigma().enable();

    item_A0Re().setState(settings().biotSavart);
    item_A0Im().setState(settings().biotSavart);

    item_ARe().enable();
    item_AIm().enable();

    item_VReGrad().enable();
    item_VImGrad().enable();

    item_BRe().enable();
    item_BIm().enable();

    item_FL().enable();
    item_pB().enable();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::eddyCurrentApp::Manager::Regions::create(const word& ccase) const
{
    region_DEFAULT().enable();
    region_CONDUCTOR().enable();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eddyCurrentApp::Manager::Manager
(
    const argList& args,
    Time& time,
    regionFvMesh& mesh,
    bool master,
    const word& name
)
:
    solverManager<regionFvMesh>
    (
        args, time, mesh, master, name
    )
{
    if (master)
    {
        // Init regionMesh
        if (mesh.initialized())
        {
            FatalErrorIn
            (
                "Foam::eddyCurrentApp::Manager::Manager(...) : "
            )
                << "Region mesh is already initialized. This Manager "
                << "needs to initialize the region mesh by itself. "
                << "Construct the regionMesh with 'init=false'!"
                << abort(FatalError);
        }
        else
        {
            HashTable<label> regionNameHashTable;

            regionNameHashTable.insert
            (
                polyMesh::defaultRegion,
                Region::DEFAULT
            );

            regionNameHashTable.insert
            (
                word(this->regionsDict().lookup("conductor")),
                Region::CONDUCTOR
            );

            wordList regionNames =
                wordList
                (
                    regionNameHashTable.size(),
                    word("MISSING_REGION_NAME")
                );

            forAllConstIter(HashTable<label>, regionNameHashTable, iter)
            {
                regionNames[iter()] = iter.key();
            }

            mesh.init(regionNames);

            this->messages().newLine();
        }
    }

    this->read();
    this->init();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eddyCurrentApp::Manager::~Manager()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::eddyCurrentApp::Manager::next() const
{}


void Foam::eddyCurrentApp::Manager::write() const
{}


void Foam::eddyCurrentApp::Manager::finalize() const
{}


// ************************************************************************* //

