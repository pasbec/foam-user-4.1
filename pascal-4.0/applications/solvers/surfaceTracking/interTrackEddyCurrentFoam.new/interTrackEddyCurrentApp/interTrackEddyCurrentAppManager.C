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

#include "interTrackEddyCurrentAppManager.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::interTrackEddyCurrentApp::Manager, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::interTrackEddyCurrentApp::Manager::Settings::read() const
{
    Foam::interTrackEddyCurrentApp::Manager::debug =
        dict().lookupOrDefault
        (
            "debug",
            Foam::interTrackEddyCurrentApp::Manager::debug()
        );

    Foam::interTrackEddyCurrentApp::Control::debug =
        dict().lookupOrDefault
        (
            "debug",
            Foam::interTrackEddyCurrentApp::Control::debug()
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// TODO: Move to conductor region
void Foam::interTrackEddyCurrentApp::Manager::Storage::Item_emPrevC::create
(const word& ccase) const
{
    set
    (
        new volVectorField
        (
            IOobject
            (
                "emPrevC",
                time().timeName(),
                mesh()[Region::CONDUCTOR],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh()[Region::CONDUCTOR],
            dimensionedVector
            (
                word(),
                dimLength,
                vector::zero
            ),
            zeroGradientFvPatchVectorField::typeName
        )
    );
}


// TODO: Move to conductor region
void Foam::interTrackEddyCurrentApp::Manager::Storage::Item_emRelDeltaA::create
(const word& init) const
{
    set
    (
        new volScalarField
        (
            IOobject
            (
                "emRelDeltaA",
                time().timeName(),
                mesh()[Region::CONDUCTOR],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh()[Region::CONDUCTOR],
            dimensionedScalar
            (
                word(),
                dimless,
                0.0
            ),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

void Foam::interTrackEddyCurrentApp::Manager::Storage::create
(
    const word& ccase
) const
{
    item_emPrevC().enable();
    item_emRelDeltaA().enable();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::interTrackEddyCurrentApp::Manager::Regions::create
(
    const word& ccase
) const
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interTrackEddyCurrentApp::Manager::Manager
(
    const argList& args,
    Time& time,
    regionDynamicFvMesh& mesh,
    bool master,
    const word& name
)
:
    solverManager<regionDynamicFvMesh>
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

            regionNameHashTable.insert
            (
                word(this->regionsDict().lookup("fluid")),
                Region::FLUID
            );

            regionNameHashTable.insert
            (
                word(this->regionsDict().lookup("buffer")),
                Region::BUFFER
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

Foam::interTrackEddyCurrentApp::Manager::~Manager()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::interTrackApp::Manager&
Foam::interTrackEddyCurrentApp::Manager::interTrackAppManager() const
{
    if (interTrackAppManager_.empty())
    {
        interTrackAppManager_.set
        (
            new interTrackApp::Manager
            (
                this->args(),
                this->time(),
                this->mesh()[Region::FLUID],
                false
            )
        );
    }

    return interTrackAppManager_();
}

Foam::eddyCurrentApp::Manager&
Foam::interTrackEddyCurrentApp::Manager::eddyCurrentAppManager() const
{
    if (eddyCurrentAppManager_.empty())
    {
        eddyCurrentAppManager_.set
        (
            new eddyCurrentApp::Manager
            (
                this->args(),
                this->time(),
                this->mesh(),
                false
            )
        );
    }

    return eddyCurrentAppManager_();
}


bool Foam::interTrackEddyCurrentApp::Manager::setCoNum
(
    scalar& CourantNumber
) const
{
    return interTrackAppManager().setCoNum(CourantNumber);
}


void Foam::interTrackEddyCurrentApp::Manager::init() const
{
    interTrackAppManager().init();
    eddyCurrentAppManager().init();
    storage().checkInit();
    regions().checkInit();
}


void Foam::interTrackEddyCurrentApp::Manager::read() const
{
    interTrackAppManager().read();
    eddyCurrentAppManager().read();
    settings().checkRead();
    regions().checkRead();
}


void Foam::interTrackEddyCurrentApp::Manager::next() const
{
    interTrackAppManager().next();
    eddyCurrentAppManager().next();
}


void Foam::interTrackEddyCurrentApp::Manager::write() const
{
    interTrackAppManager().write();
    eddyCurrentAppManager().write();
}


void Foam::interTrackEddyCurrentApp::Manager::finalize() const
{
    interTrackAppManager().finalize();
    eddyCurrentAppManager().finalize();
}


// ************************************************************************* //

