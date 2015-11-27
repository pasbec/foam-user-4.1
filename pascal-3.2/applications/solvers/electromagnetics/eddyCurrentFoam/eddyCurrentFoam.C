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

// TODO: Derived gradient boundary condition for VRe/VIm in conductor region

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createRegionMesh.H"

#   include "createRegionFields.H"
#   include "createBaseFields.H"
#   include "createConductorFields.H"

#ifndef namespaceFoam
#define namespaceFoam
    using namespace Foam;
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    while (runTime.run())
    {
        runTime++;

        Info << "Time = " << runTime.value() << nl << endl;

        // AV iterations
        while (control.loop())
        {
            // Init gradient of V in first iteration
            if (control.firstIter())
            {
                updateGradientV(control.conductor());
            }

            // Solve for A in base region
            {
                setRegionScope(control.base());

                // Relax gradient of V
                control.relax(VReGrad);
                control.relax(VImGrad);

#               include "AEqn.H"
            }

            // Update A and sigma
            updateSigmaA(control.conductor());

            // Solve for V in conductor region
            {
                setRegionScope(control.conductor());

#               include "VEqn.H"
            }

            // Update gradient of V
            updateGradientV(control.conductor());

            control.subWrite();
        }

        // Derived fields in base region
        {
            setRegionScope(control.base());

            // Magnetic field density
            BRe == fvc::curl(ARe);
            BIm == fvc::curl(AIm);

            // Eddy current density
            jRe ==   control.omega() * sigma * AIm - sigma * VReGrad;
            jIm == - control.omega() * sigma * ARe - sigma * VImGrad;

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
