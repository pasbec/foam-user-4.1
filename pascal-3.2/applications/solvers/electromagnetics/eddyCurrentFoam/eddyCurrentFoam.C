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
#include "eddyCurrentControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"

#   include "createRegionMesh.H"
#   include "createRegionControls.H"
#   include "createRegionFields.H"

#   include "createBaseFields.H"
#   include "createConductorFields.H"

    using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    while (runTime.run())
    {
        runTime++;

        Info << "Time = " << runTime.value() << nl << endl;

        // AV iterations
        while (eddy.run())
        {
            // Init gradient of V in first iteration
            if (eddy.firstIter())
            {
                updateGradientV(eddy.conductor());
            }

            // Solve for A in base region
            {
                setRegion(eddy.base());

                // Relax gradient of V
                eddy.relax(VReGrad);
                eddy.relax(VImGrad);

#               include "AEqn.H"
            }

            // Update A and sigma
            updateSigmaA(eddy.conductor());

            // Solve for V in conductor region
            {
                setRegion(eddy.conductor());

#               include "VEqn.H"
            }

            // Update gradient of V
            updateGradientV(eddy.conductor());

            eddy.subWrite();
        }

        // Derived fields in base region
        {
            setRegion(eddy.base());

            // Magnetic field density
            BRe == fvc::curl(ARe);
            BIm == fvc::curl(AIm);

            // Eddy current density
            jRe ==   eddy.omega() * sigma * AIm - sigma * VReGrad;
            jIm == - eddy.omega() * sigma * ARe - sigma * VImGrad;

            // Time-averaged Lorentz-force
            FL == 0.5 * ( (jRe ^ BRe) + (jIm ^ BIm) );

            // Time-averaged magnetic pressure
            pB == 0.5 * physicalConstant::rMu0
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
