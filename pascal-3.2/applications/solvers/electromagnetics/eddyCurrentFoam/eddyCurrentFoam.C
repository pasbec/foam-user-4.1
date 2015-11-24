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
#include "physicalConstants.H"
#include "regionModelling.H"
#include "eddyCurrentControl.H"

#include "fixedGradientFvPatchFields.H"
#include "tangentialMagneticFvPatchFields.H" // TODO: add A in name

// TODO
// #include "deflatedPCG.H"
// #include "deflatedICCG.H"

// TODO: Finish parallel face mapping
// TODO: Find bug in mapCopyInternal() member function of regionVolFields
// TODO: Finish template specialization for regionVolFields mapping

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

    // Init gradient of V
    {
        // Calculate gradient of V in conductor region
        {
            setRegionReferences(eddy.conductor());

            VReGrad = fvc::grad(VRe);
            VImGrad = fvc::grad(VIm);
        }

        // Reverse map gradient of V from conductor to base region
        VReGrad_.rmap(eddy.conductor());
        VImGrad_.rmap(eddy.conductor());
    }

    // TODO [Low]: Realize AV-loops as time loops
    while (runTime.run())
    {
        runTime++;

        Info << "Time = " << runTime.value() << nl << endl;

        // AV iterations
        while (eddy.run())
        {
            // Solve for A in base region
            {
                setRegionReferences(eddy.base());

                eddy.relax(VReGrad);
                eddy.relax(VImGrad);

#               include "AEqn.H"
            }

            // Map internal field from base to conductor region
            // and interpolate/extrapolate boundary values
            {
                ARe_.mapInterpolate(eddy.conductor());
                AIm_.mapInterpolate(eddy.conductor());
                sigma_.mapExtrapolate(eddy.conductor());
            }

            // Solve for V in conductor region
            {
                setRegionReferences(eddy.conductor());

#               include "VEqn.H"
            }

            // Gradient of V
            {
                // Calculate gradient of V in conductor region
                {
                    setRegionReferences(eddy.conductor());

                    VReGrad = fvc::grad(VRe);
                    VImGrad = fvc::grad(VIm);
                }

                // Reverse map gradient of V from conductor to base region
                VReGrad_.rmap(eddy.conductor());
                VImGrad_.rmap(eddy.conductor());

                // Store old gradient of V in base region
                {
                    setRegionReferences(eddy.base());

                    VReGrad.storePrevIter();
                    VImGrad.storePrevIter();
                }
            }

            eddy.subWrite();
        }

        // Derived fields in base region
        {
            setRegionReferences(eddy.base());

            // Magnetic field density
            BRe == fvc::curl(ARe);
            BIm == fvc::curl(AIm);

            // Eddy current density
            jRe ==   eddy.omega() * sigma * AIm - sigma * VReGrad;
            jIm == - eddy.omega() * sigma * ARe - sigma * VImGrad;
        }

        // Time-averaged fields in base region
        {
            setRegionReferences(eddy.base());

            FL == 0.5 * ( (jRe ^ BRe) + (jIm ^ BIm) );

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
