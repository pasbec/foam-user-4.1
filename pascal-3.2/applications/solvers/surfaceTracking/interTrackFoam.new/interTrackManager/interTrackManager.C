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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(interTrackManager, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void interTrackManager::Settings::read() const
{
    debug =
        dict().lookupOrDefault
        (
            "debug", Switch(interTrackManager::debug)
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void interTrackManager::Storage::init() const
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool interTrackManager::setCoNum(scalar& CourantNumber) const
{
    CourantNumber = 0.0;

    const Time& runTime = time();
    const fvMesh& mesh = this->mesh();

    DefaultRegion::Storage& storage = regions().defaultRegion().storage();

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


void interTrackManager::next() const
{
    DefaultRegion::Storage& storage = regions().defaultRegion().storage();

// TODO: Use const
    trackedSurface& interface = storage.interface();

        vector totalSurfaceTensionForce =
             interface.totalSurfaceTensionForce();

        vector totalForce =
            interface.totalViscousForce()
          + interface.totalPressureForce();

        Info << "Surface tension force: "
            << totalSurfaceTensionForce << endl;
        Info << "Surface total force: "
            << totalForce << endl;

// TODO TEST: Sub-mesh
        const areaScalarField& K = interface.curvature();
//        const areaScalarField& K = interface.aMesh().faceCurvatures();

        Info << "Surface curvature min: " << gMin(K)
            << " max: " << gMax(K)
            << " average: " << gAverage(K)
            << endl;
}


void interTrackManager::write() const
{
    DefaultRegion::Storage& storage = regions().defaultRegion().storage();

// TODO: Use const
    trackedSurface& interface = storage.interface();

    if (debug)
    {
        interface.writeVTK();
        interface.writeA();
        interface.writeVolA();
    }
}


void interTrackManager::finalize() const
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

interTrackManager::interTrackManager
(
    const argList& args,
    Time& time,
    dynamicFvMesh& mesh,
    const word& name,
    const bool& master
)
:
    solverManager<dynamicFvMesh>
    (
        args, time, mesh, name, master
    )
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

