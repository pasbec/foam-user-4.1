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

#include "cellSet.H"
#include "faceSet.H"

#include "fvMeshSubset.H"

#include "electricPotentialGrad.H"
#include "electricPotentialLaplacian.H"

#include "zeroGradientFvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "lookupFixedGradientFvPatchFields.H"
#include "tangentialMagneticFvPatchFields.H"

#include "deflatedPCG.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "addRegionOption.H"

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createNamedMesh.H"

    using namespace Foam;

#   include "createFields.H"
#   include "createTopo.H"

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

        // Solve AV system
        {
            // Read controls
#           include "readAVControls.H"

            // Solve for A and V
#           include "solveAV.H"
        }

        // Derived fields
        {
            // Magnetic field density
            BRe == fvc::curl(ARe);
            BIm == fvc::curl(AIm);

            // Eddy current density
            jRe ==   alpha * AIm - sigma * VReGrad;
            jIm == - alpha * ARe - sigma * VImGrad;
        }

        // Time-averaged fields
        {
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
