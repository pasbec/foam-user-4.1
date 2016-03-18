/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description

\*---------------------------------------------------------------------------*/

#include "faSubMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void faSubMesh::makeFaceCurvatures() const
{
    if (debug)
    {
        Info<< "faSubMesh::makeFaceCurvatures() : "
            << "making surface normal derivative of normal velocity"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (faceCurvaturesPtr_)
    {
        FatalErrorIn("faSubMesh::makeFaceCurvatures()")
            << "face curvatures field already exists"
                << abort(FatalError);
    }

    faceCurvaturesPtr_ =
        new areaScalarField
        (
            IOobject
            (
                "faceCurvatures",
                basePolyMesh().pointsInstance(),
                faMesh::meshSubDir,
                basePolyMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            baseAreaMesh(),
            dimensionedScalar
            (
                word(),
                dimless/dimLength,
                0
            )
        );

    areaScalarField& K = *faceCurvaturesPtr_;

    scalarField& KIn = K.internalField();

    const areaScalarField& subK = subAreaMesh().faceCurvatures();

    const pointField& subPoints = subAreaMesh().points();
    const vectorField& subFaceCentres = subAreaMesh().areaCentres();
    const labelListList& subPointFaces = subAreaMesh().patch().pointFaces();

    forAll (KIn, facei)
    {
        label subPointI = baseAreaMesh().nPoints() + facei;

// TODO: This might fail for some cases!
// WARNING
//
// subPointI is a GLOBAL point index, but
// subPointFaces takes LOCAL point indices from
// this patch!!!
// Use baseAreaMesh().patch().meshPointMap() which is SLOW!!!!
// Or find some other way..

        const labelList& subFacesI = subPointFaces[subPointI];

        scalar wSum = 0.0;

        forAll (subFacesI, subFacei)
        {
            label subFaceI = subFacesI[subFacei];

            scalar wI = 1.0
                / (mag(subFaceCentres[subFaceI]-subPoints[subPointI]) + SMALL);

            wSum += wI;

            KIn[facei] += wI * subK[subFaceI];

        }

        KIn[facei] /= wSum + SMALL;
    }
// Pout << "K.boundaryField() = " << K.boundaryField() << endl;
// Pout << "subK.boundaryField() = " << subK.boundaryField() << endl;
// TODO: This is totally wrong!
//     K.boundaryField() == subK.boundaryField();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const areaScalarField& faSubMesh::faceCurvatures() const
{
    if (!faceCurvaturesPtr_)
    {
        makeFaceCurvatures();
    }

    return *faceCurvaturesPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
