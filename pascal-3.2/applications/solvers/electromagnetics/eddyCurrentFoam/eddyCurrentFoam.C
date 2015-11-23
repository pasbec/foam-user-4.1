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

#include "wordIOList.H"

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
#   include "createRegionFields.H"

    // Eddy current controls
    eddyCurrentControl eddy(control_);

#   include "createBaseFields.H"
#   include "readBaseControls.H"
#   include "initBaseControls.H"

#   include "createConductorFields.H"
#   include "readConductorControls.H"
#   include "initConductorControls.H"

    using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Initial gradient calculation of V in conductor region
    Info << "Init gradient of V" << nl << endl;
    {
        setRegionReferences(eddy.conductorRegion());

        VReGrad = fvc::grad(VRe);
        VImGrad = fvc::grad(VIm);
    }

    // Initial reverse map gradient of V from conductor to base region
    VReGrad_.rmap(eddy.conductorRegion());
    VImGrad_.rmap(eddy.conductorRegion());

    // TODO [Low]: Realize AV-loops as time loops
    while (runTime.run())
    {
        runTime++;

        Info << "Time = " << runTime.value() << nl << endl;

        // Copy local solver dictionaries
        dictionary Adict = eddy.Adict();
        dictionary Vdict = eddy.Vdict();

        // Reset current tolerance
        curTol = 1.0;

        // AV iterations
        for (int iter = 0; iter < eddy.maxIter(); iter++)
        {
            if ((iter >= eddy.minIter()) && (res < eddy.tol())) break;

            if (res < curTol)
            {
                curTol /= 10.0;

                Info << nl << "Tolerance = " << curTol
                    << " / " << eddy.tol() << endl << nl;

                Adict.set<scalar>("tolerance", curTol);
                Vdict.set<scalar>("tolerance", curTol);
            }

            // Solve for A in base region
            {
                setRegionReferences(eddy.baseRegion());

#               include "AEqn.H"
            }

            // Map internal field from base to conductor region
            // and interpolate/extrapolate boundary values
            ARe_.mapInterpolate(eddy.conductorRegion());
            AIm_.mapInterpolate(eddy.conductorRegion());
            sigma_.mapExtrapolate(eddy.conductorRegion());

            // Solve for V in conductor region
            {
                setRegionReferences(eddy.conductorRegion());

#               include "VEqn.H"
            }

            // Calculate gradient of V in conductor region
            {
                setRegionReferences(eddy.conductorRegion());

                VReGrad = fvc::grad(VRe);
                VImGrad = fvc::grad(VIm);
            }

            // Reverse map gradient of V from conductor to base region
            VReGrad_.rmap(eddy.conductorRegion());
            VImGrad_.rmap(eddy.conductorRegion());

            // Update residual
            res = max(resA, resV);
            reduce(res, maxOp<scalar>());

            if (res < curTol)
            {
                runTime.writeNow();
            }
        }

        // Derived fields in base region
        {
            setRegionReferences(eddy.baseRegion());

            // Magnetic field density
            BRe == fvc::curl(ARe);
            BIm == fvc::curl(AIm);

            // Eddy current density
            jRe ==   eddy.omega() * sigma * AIm - sigma * VReGrad;
            jIm == - eddy.omega() * sigma * ARe - sigma * VImGrad;
        }

        // Time-averaged fields in base region
        {
            setRegionReferences(eddy.baseRegion());

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
