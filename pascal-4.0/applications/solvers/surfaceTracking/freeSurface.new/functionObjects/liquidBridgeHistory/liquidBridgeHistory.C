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

#include "liquidBridgeHistory.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "boundBox.H"
#include "faMesh.H"
#include "freeSurface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(liquidBridgeHistory, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        liquidBridgeHistory,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liquidBridgeHistory::liquidBridgeHistory
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
//     historyFilePtr_(NULL),
    freeSurfacePatchID_(-1)
{
    if (Pstream::parRun())
    {
        FatalErrorIn("torusFilmHistory::torusFilmHistory(...)")
            << "torusFilmHistory objec function "
                << "is not implemented for parallel run"
                << abort(FatalError);
    }

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

//     if (Pstream::master())
//     {
//         if (!time_.processorCase())
//         {
//             mkDir
//             ( 
//                 time_.path()
//                /"history"
//                /time_.timeName()
//             );
        
//             historyFilePtr_ = 
//                 new OFstream
//                 (
//                     time_.path()
//                    /"history"
//                    /time_.timeName()
//                    /"history.dat"
//                 );
//         }
//         else
//         {
//             mkDir
//             ( 
//                 time_.path()/".."/"history"
//                /time_.timeName()
//             );
        
//             historyFilePtr_ = 
//                 new OFstream
//                 (
//                     time_.path()/".."
//                    /"history"
//                    /time_.timeName()
//                    /"history.dat"
//                 );
//         }

//         (*historyFilePtr_) 
//             << "Time" << tab 
//                 << "hLeft" << tab
//                 << "hRight" << endl;
//     }

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    freeSurfacePatchID_ = 
        mesh.boundaryMesh().findPatchID("freeSurface");

    if (freeSurfacePatchID_ == -1)
    {
        FatalErrorIn("torusFilmHistory::torusFilmHistory(...)")
            << "Problem with finding freeSurface patch"
                << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::liquidBridgeHistory::start()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const freeSurface& fs = 
        mesh.lookupObject<freeSurface>("freeSurfaceProperties");

    const vectorField& Cf = 
        fs.aMesh().areaCentres().internalField();

    const scalarField& Kf = 
        fs.aMesh().faceCurvatures().internalField();

    const scalarField& Tf = 
        fs.temperature().internalField();

    scalarField gradTf =
        (
            vector(0,1,0)
          & fac::grad(fs.temperature())().internalField()
        );

    if (Pstream::master())
    {
        {
            OFstream file
            (
                time_.timePath()/"kFs.dat"
            );

            file.precision(12);

            forAll(Cf, faceI)
            {
                file << Cf[faceI].y() << tab << Kf[faceI] << tab
                    << sqrt(sqr(Cf[faceI].x())+sqr(Cf[faceI].z())) << tab
                    << Tf[faceI] << tab << gradTf[faceI] << endl;
            }
        }

        return true;
    }

    return false;
}


bool Foam::liquidBridgeHistory::execute()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const freeSurface& fs = 
        mesh.lookupObject<freeSurface>("freeSurfaceProperties");

    const vectorField& Cf = 
        fs.aMesh().areaCentres().internalField();

    const scalarField& Kf = 
        fs.aMesh().faceCurvatures().internalField();
    
    const scalarField& Tf = 
        fs.temperature().internalField();

    scalarField gradTf =
        (
            vector(0,1,0)
          & fac::grad(fs.temperature())().internalField()
        );

    if (Pstream::master())
    {
        if (time_.outputTime())
        {
            OFstream file
            (
                time_.timePath()/"kFs.dat"
            );

            file.precision(12);

            forAll(Cf, faceI)
            {
                file << Cf[faceI].y() << tab << Kf[faceI] << tab
                    << sqrt(sqr(Cf[faceI].x())+sqr(Cf[faceI].z())) << tab
                    << Tf[faceI] << tab << gradTf[faceI] << endl;
            }
        }

        return true;
    }

    return false;
}


bool Foam::liquidBridgeHistory::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}

// ************************************************************************* //
