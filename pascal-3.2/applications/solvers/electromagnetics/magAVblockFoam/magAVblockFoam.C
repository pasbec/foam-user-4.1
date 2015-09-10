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

Application
    ...

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvBlockMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    using namespace Foam;

#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Constants, Frequency and Alpha

    dimensionedScalar mu0
    (
        "mu0",
        dimensionSet(1, 1, -2, 0, 0, -2, 0),
        4.0 * mathematicalConstant::pi * 1e-07
    );

    dimensionedScalar rMu0
    (
        "rMu0",
        1.0/mu0
    );

    dimensionedScalar omega
    (
        "omega",
        2.0 * mathematicalConstant::pi * frequency
    );

    volScalarField alpha
    (
        "alpha",
        omega * sigma
    );

    while (runTime.run())
    {
        runTime++;

        Info << "Time = " << runTime.value() << nl << endl;

        // Re-read controls
#       include "readAVControls.H"

        // Solve AV-system
#       include "solveAV.H"

        jRe ==   alpha * AIm - sigma * fvc::grad(VRe);
        jIm == - alpha * ARe - sigma * fvc::grad(VIm);

        BRe == fvc::curl(ARe);
        BIm == fvc::curl(AIm);

        FL == 0.5 * ( (jRe ^ BRe) + (jIm ^ BIm) );

        pB == 0.5 * rMu0
            * 0.5 * ( (BRe & BRe) + (BIm & BIm) );

        // ==================================================================//
        // Debug stuff
        // ==================================================================//

        gradAlpha.write();

        surfaceScalarField AReFlux ("AReFlux", fvc::interpolate(ARe) & mesh.Sf());
        surfaceScalarField AImFlux ("AImFlux", fvc::interpolate(AIm) & mesh.Sf());

        volScalarField AReMag ("AReMag", mag(ARe)); AReMag.write();
        volScalarField AImMag ("AImMag", mag(AIm)); AImMag.write();

        volVectorField AReRes
        (
            "AReRes",
            rMu0 * fvc::laplacian(ARe)
          + alpha * AIm
          - sigma * fvc::grad(VRe)
          + jsRe
        ); AReRes.write();

        volVectorField AImRes
        (
            "AImRes",
            rMu0 * fvc::laplacian(AIm)
          - alpha * ARe
          - sigma * fvc::grad(VIm)
          + jsIm
        ); AImRes.write();

        volScalarField VReRes
        (
            "VReRes",
            fvc::laplacian(sigmaf, VRe)
          - (gradAlpha&AIm)
        ); VReRes.write();

        volScalarField VImRes
        (
            "VImRes",
            fvc::laplacian(sigmaf, VIm)
          + (gradAlpha&ARe)
        ); VImRes.write();

        volScalarField BReMag ("BReMag", mag(BRe)); BReMag.write();
        volScalarField BImMag ("BImMag", mag(BIm)); BImMag.write();

        surfaceScalarField jReFlux ("jReFlux", fvc::interpolate(jRe) & mesh.Sf());
        surfaceScalarField jImFlux ("jImFlux", fvc::interpolate(jIm) & mesh.Sf());

        volScalarField jReMag ("jReMag", mag(jRe)); jReMag.write();
        volScalarField jImMag ("jImMag", mag(jIm)); jImMag.write();

        volScalarField jsReMag ("jsReMag", mag(jsRe)); jsReMag.write();
        volScalarField jsImMag ("jsImMag", mag(jsIm)); jsImMag.write();

        // ==================================================================//
        // Write
        // ==================================================================//

        runTime.write();

        Info << "ExecutionTime = "
            << scalar(runTime.elapsedCpuTime()) << " s"
            << nl << endl << endl;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "End" << nl << endl;

    return(0);
}

// ************************************************************************* //
