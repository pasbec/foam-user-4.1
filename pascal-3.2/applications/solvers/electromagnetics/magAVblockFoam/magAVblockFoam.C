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
#include "faceSet.H"

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

    dimensionedScalar rMu0(
        "rMu0",
        1.0/mu0
    );

    dimensionedScalar omega
    (
        "omega",
        2.0 * mathematicalConstant::pi * frequency
    );

    dimensionedScalar rOmega(
        "rOmega",
        1.0/omega
    );

    volScalarField alpha(
        "alpha",
        omega * sigma
    );

    while (runTime.run())
    {
        runTime++;

        Info << "Time = " << runTime.value() << nl << endl;

        // ==================================================================//
        // Solve quasi-static Maxwell-Equations for low Rm
        // ==================================================================//

        // Interpolate sigma to face centers
        surfaceScalarField sigmaf(
            "sigmaf",
            fvc::interpolate(sigma,"interpolate(sigma)")
        );

        // Get alpha on the face centers
        surfaceScalarField alphaf(
            "alphaf",
            omega*fvc::interpolate(sigma,"interpolate(alpha)")

            );

        // Prepare block system
        fvBlockMatrix<vector10> AVEqn(AV);

        // Assembleand insert equations for A
#       include "AEqn.H"

        AVEqn.insertEquation(0, AReEqn);
        AVEqn.insertEquation(5, AImEqn);

        // Assembleand insert equations for divA
#       include "divAEqn.H"

        AVEqn.insertEquation(3, divAReEqn);
        AVEqn.insertEquation(8, divAImEqn);

        // Assemble and insert equations for V
#       include "VEqn.H"

        AVEqn.insertEquation(4, VReEqn);
        AVEqn.insertEquation(9, VImEqn);

        // Assemble and insert coupling terms
#       include "couplingTerms.H"

        //- Block coupled solver call
        AVEqn.solve();

        // Retrieve solution
        AVEqn.retrieveSolution(0, ARe.internalField());
        AVEqn.retrieveSolution(5, AIm.internalField());

        ARe.correctBoundaryConditions();
        AIm.correctBoundaryConditions();

        AVEqn.retrieveSolution(3, divARe.internalField());
        AVEqn.retrieveSolution(8, divAIm.internalField());

        ARe.correctBoundaryConditions();
        AIm.correctBoundaryConditions();

        AVEqn.retrieveSolution(4, VRe.internalField());
        AVEqn.retrieveSolution(9, VIm.internalField());

        VRe.correctBoundaryConditions();
        VIm.correctBoundaryConditions();

        // ==================================================================//
        // Calculate derived fields
        // ==================================================================//

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

        surfaceScalarField AReFlux ("AReFlux", fvc::interpolate(ARe) & mesh.Sf());
        surfaceScalarField AImFlux ("AImFlux", fvc::interpolate(AIm) & mesh.Sf());

        volScalarField AReMag ("AReMag", mag(ARe)); AReMag.write();
        volScalarField AImMag ("AImMag", mag(AIm)); AImMag.write();

        volScalarField BReMag ("BReMag", mag(BRe)); BReMag.write();
        volScalarField BImMag ("BImMag", mag(BIm)); BImMag.write();

        surfaceScalarField jReFlux ("AReFlux", fvc::interpolate(jRe) & mesh.Sf());
        surfaceScalarField jImFlux ("AImFlux", fvc::interpolate(jIm) & mesh.Sf());

        volScalarField jReMag ("jReMag", mag(jRe)); jReMag.write();
        volScalarField jImMag ("jImMag", mag(jIm)); jImMag.write();

        volScalarField jReDiv ("jReDiv", fvc::div(jRe)); jReDiv.write();
        volScalarField jImDiv ("jImDiv", fvc::div(jIm)); jImDiv.write();

        volScalarField jsReMag ("jsReMag", mag(jsRe)); jsReMag.write();
        volScalarField jsImMag ("jsImMag", mag(jsIm)); jsImMag.write();

        faceSet conductorFaces(mesh, "faceSet_geometry_conductor");

        scalarField conductorSigmaf(conductorFaces.size(),0.0);
        scalarField conductorAlphaAReFlux(conductorFaces.size(),0.0);
        scalarField conductorAlphaAImFlux(conductorFaces.size(),0.0);
        scalarField conductorJReFlux(conductorFaces.size(),0.0);
        scalarField conductorJImFlux(conductorFaces.size(),0.0);

        forAll(conductorFaces, faceI)
        {
            label faceLabel = conductorFaces[faceI];

            conductorSigmaf[faceI] = sigmaf[faceLabel];

            conductorAlphaAReFlux[faceI] = omega.value() * sigmaf[faceI] * AReFlux[faceI];
            conductorAlphaAImFlux[faceI] = omega.value() * sigmaf[faceI] * AImFlux[faceI];

            conductorJReFlux[faceI] = jReFlux[faceLabel];
            conductorJImFlux[faceI] = jImFlux[faceLabel];
        }

        // ==================================================================//
        // Debug output
        // ==================================================================//

        Info << "divAReEqn.hasDiag()" << " = " << divAReEqn.hasDiag() << endl;
        Info << "divAReEqn.hasLower()" << " = " << divAReEqn.hasLower() << endl;
        Info << "divAReEqn.hasUpper()" << " = " << divAReEqn.hasUpper() << endl;
        Info << "divAReEqn.diagonal()" << " = " << divAReEqn.diagonal() << endl;

        Info << "divAImEqn.hasDiag()" << " = " << divAImEqn.hasDiag() << endl;
        Info << "divAImEqn.hasLower()" << " = " << divAImEqn.hasLower() << endl;
        Info << "divAImEqn.hasUpper()" << " = " << divAImEqn.hasUpper() << endl;
        Info << "divAImEqn.diagonal()" << " = " << divAImEqn.diagonal() << endl;

        Info << "max(mag(jReDiv)) = " << max(mag(jReDiv)) << endl;
        Info << "max(mag(jImDiv)) = " << max(mag(jImDiv)) << endl;

        Info << "gMax(conductorSigmaf) = " << gMax(conductorSigmaf) << endl;
        Info << "gMax(mag(conductorAlphaAReFlux)) = " << gMax(mag(conductorAlphaAReFlux)) << endl;
        Info << "gMax(mag(conductorAlphaAImFlux)) = " << gMax(mag(conductorAlphaAImFlux)) << endl;
        Info << "gMax(mag(conductorJReFlux)) = " << gMax(mag(conductorJReFlux)) << endl;
        Info << "gMax(mag(conductorJImFlux)) = " << gMax(mag(conductorJImFlux)) << endl;

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
