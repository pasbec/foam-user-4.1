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
#include "harmonic.H"

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

    dimensionedScalar rMu0 = 1.0/mu0;

    dimensionedScalar omega
    (
        "omega",
        2.0 * mathematicalConstant::pi * frequency
    );
    
    dimensionedScalar rOmega = 1.0/omega;

    volScalarField alpha = omega * sigma;
        
    while (runTime.run())
    {
        runTime++;

        Info << "Time = " << runTime.value() << nl << endl;
        
//         // Harmonic interpolation schemes
//         harmonic<scalar> harmonicScalarScheme(mesh);
//         harmonic<vector> harmonicVectorScheme(mesh);
    
        // Interpolate sigma to face centers
        surfaceScalarField sigmaf("sigmaf", fvc::interpolate(sigma,"interpolate(sigma)"));
    
        // Get alpha on the face centers
        surfaceScalarField alphaf("alphaf", fvc::interpolate(omega*sigma,"interpolate(alpha)"));
        
        // ==================================================================//
        // Solve quasi-static Maxwell-Equations for low Rm
        // ==================================================================//
        
        // Prepare block system
        fvBlockMatrix<vector8> AVEqn(AV);
        
        // Assemble equations for A
#       include "AEqn.H"

        // Assemble and insert pressure equation
#       include "VEqn.H"

        // Insert equations into block Matrix
        AVEqn.insertEquation(0, AReEqn);
        AVEqn.insertEquation(4, AImEqn);
        
        AVEqn.insertEquation(3, VReEqn);
        AVEqn.insertEquation(7, VImEqn);

        // Assemble and insert coupling terms
#       include "couplingTerms.H"

        //- Block coupled solver call
        AVEqn.solve();

        // Retrieve solution
        AVEqn.retrieveSolution(0, ARe.internalField());
        AVEqn.retrieveSolution(4, AIm.internalField());

        ARe.correctBoundaryConditions();
        AIm.correctBoundaryConditions();
        
        AVEqn.retrieveSolution(3, VRe.internalField());
        AVEqn.retrieveSolution(7, VIm.internalField());

        VRe.correctBoundaryConditions();
        VIm.correctBoundaryConditions();

        // ==================================================================//
        // Calculate derived fields
        // ==================================================================//

        jRe ==  alpha * AIm;
        jIm == -alpha * ARe;

        BRe == fvc::curl(ARe);
        BIm == fvc::curl(AIm);

        FL == 0.5 * ( (jRe ^ BRe) + (jIm ^ BIm) );

        pB == 0.5 * rMu0
            * 0.5 * ( (BRe & BRe) + (BIm & BIm) );

        // ==================================================================//
        // Debug fields
        // ==================================================================//
        
        volScalarField AReMag ("AReMag", mag(ARe)); AReMag.write();
        volScalarField AImMag ("AImMag", mag(AIm)); AImMag.write();
        
        volScalarField BReMag ("BReMag", mag(BRe)); BReMag.write();
        volScalarField BImMag ("BImMag", mag(BIm)); BImMag.write();
        
        volScalarField jReMag ("jReMag", mag(jRe)); jReMag.write();
        volScalarField jImMag ("jImMag", mag(jIm)); jImMag.write();
        
        volScalarField jReDiv ("jReDiv", fvc::div(jRe)); jReDiv.write();
        volScalarField jImDiv ("jImDiv", fvc::div(jIm)); jImDiv.write();
        
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
