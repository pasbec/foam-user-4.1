/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*----------------------------------------------------------------------------*/

#include "dropletHistory.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "boundBox.H"
#include "trackedSurface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dropletHistory, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        dropletHistory,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dropletHistory::dropletHistory
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    historyFilePtr_(NULL),
    trackedSurfacePatchID_(-1),
    V0_(SMALL)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    if (Pstream::master())
    {
        if (!time_.processorCase())
        {
            mkDir
            (
                time_.path()
               /"history"
               /time_.timeName()
            );

            historyFilePtr_ =
                new OFstream
                (
                    time_.path()
                   /"history"
                   /time_.timeName()
                   /"history.dat"
                );
        }
        else
        {
            mkDir
            (
                time_.path()/".."/"history"
               /time_.timeName()
            );

            historyFilePtr_ =
                new OFstream
                (
                    time_.path()/".."
                   /"history"
                   /time_.timeName()
                   /"history.dat"
                );
        }

        (*historyFilePtr_)
            << "Time" << tab
                << "Cx" << tab
                << "Cy" << tab
                << "Cz" << tab
                << "Xmin" << tab
                << "Xmax" << tab
                << "Ymin" << tab
                << "Ymax" << tab
                << "Zmin" << tab
                << "Zmax" << tab
                << "V/V0" << tab
                << "Fx" << tab
                << "Fy" << tab
                << "Fz" << endl;
    }

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    trackedSurfacePatchID_ = mesh.boundaryMesh().findPatchID("trackedSurface");

    if (trackedSurfacePatchID_ == -1)
    {
        FatalErrorIn("dropHistory::dropHistory(...)")
            << "Problem with finding trackedSurface patch"
                << abort(FatalError);
    }

    V0_ = gSum(mesh.V());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::dropletHistory::start()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    vector C = gSum(mesh.C().internalField()*mesh.V())/gSum(mesh.V());

    boundBox box(mesh.C().boundaryField()[trackedSurfacePatchID_]);

    scalar V = gSum(mesh.V());

    const trackedSurface& ts =
        mesh.lookupObject<trackedSurface>("trackedSurfaceProperties");

    vector totalForce = ts.totalViscousForce() + ts.totalPressureForce();

    if (Pstream::master())
    {
        historyFilePtr_->precision(12);

        (*historyFilePtr_) << time_.value() << tab
            << C.x() << tab
            << C.y() << tab
            << C.z() << tab
            << box.min().x() << tab
            << box.max().x() << tab
            << box.min().y() << tab
            << box.max().y() << tab
            << box.min().z() << tab
            << box.max().z() << tab
            << mag(1.0 - V/V0_) << tab
            << totalForce.x() << tab
            << totalForce.y() << tab
            << totalForce.z() << endl;

        return true;
    }

    return false;
}


bool Foam::dropletHistory::execute()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    vector C = gSum(mesh.C().internalField()*mesh.V())/gSum(mesh.V());

    boundBox box(mesh.C().boundaryField()[trackedSurfacePatchID_]);

    scalar V = gSum(mesh.V());

    const trackedSurface& ts =
        mesh.lookupObject<trackedSurface>("trackedSurfaceProperties");

    vector totalForce = ts.totalViscousForce() + ts.totalPressureForce();

    if (Pstream::master())
    {
        historyFilePtr_->precision(12);

        (*historyFilePtr_) << time_.value() << tab
            << C.x() << tab
            << C.y() << tab
            << C.z() << tab
            << box.min().x() << tab
            << box.max().x() << tab
            << box.min().y() << tab
            << box.max().y() << tab
            << box.min().z() << tab
            << box.max().z() << tab
            << mag(1.0 - V/V0_) << tab
            << totalForce.x() << tab
            << totalForce.y() << tab
            << totalForce.z() << endl;

        return true;
    }

    return false;
}


bool Foam::dropletHistory::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}

// ************************************************************************* //
