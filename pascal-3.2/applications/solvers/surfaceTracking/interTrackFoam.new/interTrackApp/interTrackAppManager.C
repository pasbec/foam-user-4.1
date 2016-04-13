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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::interTrackApp::Manager, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::interTrackApp::Manager::Settings::read() const
{
    debug = dict().lookupOrDefault("debug", Switch(debug));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::interTrackApp::Manager::Storage::create(const word& ccase) const
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::interTrackApp::Manager::Regions::create(const word& ccase) const
{
    region_DEFAULT().enable();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::interTrackApp::Manager::setCoNum(scalar& CourantNumber) const
{
    CourantNumber = 0.0;

    const Time& runTime = time();
    const fvMesh& mesh = this->mesh();

    Region_DEFAULT::Storage& storage = regions().region_DEFAULT().storage();

    tmp<surfaceScalarField> tphi(storage.phi());
    surfaceScalarField& phi = tphi();

    // Convective Courant Number
    {
        if (mesh.moving())
        {
            const volVectorField& U = storage.U();

            // Make fluxes relative
            phi -= fvc::meshPhi(U);
        }

#       include "CourantNo.H"

        CourantNumber = max(CourantNumber, CoNum);
    }

    // Mesh Courant Number
    {
#       include "meshCourantNo.H"

        CourantNumber = max(CourantNumber, meshCoNum);
    }

    // Interface Courant Number
    {
// TODO: Use const
        trackedSurface& interface = storage.interface();

        scalar interfaceCoNum = interface.maxCourantNumber();

        Info << "Surface Courant Number max: "
            << interface.maxCourantNumber() << endl;

        CourantNumber = max(CourantNumber, interfaceCoNum);
    }

    tphi.clear();

    return true;
}


void Foam::interTrackApp::Manager::next() const
{
    Region_DEFAULT::Storage& storage = regions().region_DEFAULT().storage();

// TODO: Use const
    trackedSurface& interface = storage.interface();

// TODO TEST: Sub-mesh
#   include "trackedSurfaceStatistics.H"
}


void Foam::interTrackApp::Manager::write() const
{
    Region_DEFAULT::Storage& storage = regions().region_DEFAULT().storage();

// TODO: Use const
    trackedSurface& interface = storage.interface();

    if (debug)
    {
        interface.writeVTK();
        interface.writeA();
        interface.writeVolA();
    }
}


void Foam::interTrackApp::Manager::finalize() const
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interTrackApp::Manager::Manager
(
    const argList& args,
    Time& time,
    dynamicFvMesh& mesh,
    const word& name,
    const bool& master
)
:
    solverManager<dynamicFvMesh, Region::SIZE>
    (
        args, time, mesh, name, master
    )
{
    // Region name list
    regionNames_.setSize(Region::SIZE);
    regionNames_[Region::DEFAULT] = polyMesh::defaultRegion;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interTrackApp::Manager::~Manager()
{}


// ************************************************************************* //

