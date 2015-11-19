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
#include "regionModelling.H"

#include "fixedGradientFvPatchFields.H"
#include "tangentialMagneticFvPatchFields.H" // TODO: add A in name

// TODO
// #include "deflatedPCG.H"
// #include "deflatedICCG.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "addRegionOption.H"

#   include "setRootCase.H"
#   include "createTime.H"

#   include "createRegionMesh.H"
#   include "createRegionFields.H"
#   include "createRegionControl.H"

#   include "createBaseFields.H"
#   include "readBaseControls.H"
#   include "initBaseControls.H"

#   include "createConductorFields.H"
#   include "readConductorControls.H"
#   include "initConductorControls.H"

    using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Constants, Frequency and Alpha

    // Magnetic properties
    Info << "Read Magnetic properties" << nl << endl;
    IOdictionary magneticProperties
    (
        IOobject
        (
            "magneticProperties",
            runTime.constant(),
            mesh_[defaultRegionID],
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar frequency (magneticProperties.lookup("frequency"));

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

    // Calculate gradient of V
    {
        setRegionReferences(conductorRegionID);

        VReGrad = fvc::grad(VRe);
        VImGrad = fvc::grad(VIm);
    }

    // Initialize gradient of V in base region
    Info << "Init gradient of V" << nl << endl;
    VReGrad_.rmap(conductorRegionID);
    VImGrad_.rmap(conductorRegionID);

    // Initialize conductivity in conductor region
    Info << "Init sigma" << nl << endl;
    sigma_.map(conductorRegionID);

    // TODO: Realize AV-loops as time loops
    while (runTime.run())
    {
        runTime++;

        Info << "Time = " << runTime.value() << nl << endl;

        // Copy local solver dictionaries
        dictionary Adict = mesh_[defaultRegionID].solutionDict().subDict("solvers").subDict("A");
        dictionary Vdict = mesh_[conductorRegionID].solutionDict().subDict("solvers").subDict("V");

        // Reset current tolerance
        curTol = 1.0;

        // AV iterations
        for (int iter = 0; iter < maxIter; iter++)
        {
            if ((iter >= minIter) && (res < tol)) break;

            if (res < curTol)
            {
                curTol /= 10.0;

                Info << nl << "Tolerance = " << curTol
                    << " / " << tol << endl << nl;

                Adict.set<scalar>("tolerance", curTol);
                Vdict.set<scalar>("tolerance", curTol);
            }

            // Solve for A in base region
            {
                setRegionReferences(defaultRegionID);

#               include "AEqn.H"
            }

            // Map fields from base to conductor region
            ARe_.map(conductorRegionID);
            AIm_.map(conductorRegionID);
            sigma_.map(conductorRegionID);

            // Solve for V in conductor region
            {
                setRegionReferences(conductorRegionID);

#               include "VEqn.H"
            }

            // Reverse map fields from base to conductor region
            VReGrad_.rmap(conductorRegionID);
            VImGrad_.rmap(conductorRegionID);

            // Update residual
            res = max(resA, resV);
            reduce(res, maxOp<scalar>());

            if (res < curTol)
            {
                runTime.writeNow();
            }
        }

        // Derived fields
        {
            setRegionReferences(defaultRegionID);

            // Magnetic field density
            BRe == fvc::curl(ARe);
            BIm == fvc::curl(AIm);

            // Eddy current density
            jRe ==   omega * sigma * AIm - sigma * VReGrad;
            jIm == - omega * sigma * ARe - sigma * VImGrad;
        }

        // Time-averaged fields
        {
            setRegionReferences(defaultRegionID);

            FL == 0.5 * ( (jRe ^ BRe) + (jIm ^ BIm) );

            pB == 0.5 * rMu0
                * 0.5 * ( (BRe & BRe) + (BIm & BIm) );
        }

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
