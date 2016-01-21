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

#include "sloshingHistory.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "boundBox.H"
#include "faMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sloshingHistory, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        sloshingHistory,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sloshingHistory::sloshingHistory
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
    leftPointID_(-1),
    rightPointID_(-1)
{
    if (Pstream::parRun())
    {
        FatalErrorIn("sloshingHistory::sloshingHistory(...)")
            << "sloshingHistory objec function "
                << "is not implemented for parallel run"
                << abort(FatalError);
    }

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
                << "Yleft" << tab
                << "Yright" << endl;
    }

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    trackedSurfacePatchID_ =
        mesh.boundaryMesh().findPatchID("trackedSurface");

    if (trackedSurfacePatchID_ == -1)
    {
        FatalErrorIn("dropHistory::dropHistory(...)")
            << "Problem with finding trackedSurface patch"
                << abort(FatalError);
    }

    const vectorField& tsPoints =
        mesh.boundaryMesh()[trackedSurfacePatchID_].localPoints();

    scalar maxX = 0;
    scalar minX = GREAT;

    forAll(tsPoints, pointI)
    {
        if (tsPoints[pointI].x() < minX)
        {
            minX = tsPoints[pointI].x();
            leftPointID_ = pointI;
        }

        if (tsPoints[pointI].x() > maxX)
        {
            maxX = tsPoints[pointI].x();
            rightPointID_ = pointI;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sloshingHistory::start()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const vectorField& tsPoints =
        mesh.boundaryMesh()[trackedSurfacePatchID_].localPoints();

    if (Pstream::master())
    {
        historyFilePtr_->precision(12);

        (*historyFilePtr_) << time_.value() << tab
            << tsPoints[leftPointID_].y() << tab
            << tsPoints[rightPointID_].y() << endl;

        return true;
    }

    return false;
}


bool Foam::sloshingHistory::execute()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const vectorField& tsPoints =
        mesh.boundaryMesh()[trackedSurfacePatchID_].localPoints();

    if (Pstream::master() && time_.outputTime())
    {
        historyFilePtr_->precision(12);

        (*historyFilePtr_) << time_.value() << tab
            << tsPoints[leftPointID_].y() << tab
            << tsPoints[rightPointID_].y() << endl;

        return true;
    }

    return false;
}


bool Foam::sloshingHistory::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}

// ************************************************************************* //
