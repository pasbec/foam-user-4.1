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
    interTrackMHDFoam

Description

TODO [High]: Where does the very high (time-linear) volume conservation error come from if the Alberto-method is NOT used?

TODO [High]: Why is the Alberto-method much less stable in icf-small or for strong deformations - boundary correction?

TODO [High]: What about curvature at wall contact line? It seems that the contact angle is always 90°

TODO [High]: Why does tangential surface tension lead to instable calculations?

TODO [High]: What about tangential surface tension at wall contact line and with presence of FL

TODO [High]: What happens on wedge meshes and tetFemMotionSolver?

TODO [High]: What happens on annulus wedge meshes (with symmetry at annulus inner side) and pseudoSolid tetFemMotionSolver?

TODO [Low]: Add a description

\*---------------------------------------------------------------------------*/

#include <limits>

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "transportModel.H"
#include "twoPhaseMixture.H"
#include "turbulenceModel.H"
#include "trackedSurface.H"

#include "fixedGradientFvPatchFields.H"

#include "blockLduSolvers.H"
#include "VectorNFieldTypes.H"
#include "volVectorNFields.H"
#include "blockVectorNMatrices.H"

#include "blockMatrixTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createDynamicFvMesh.H"
#   include "createFields.H"
#   include "createBlockMatrixField.H"

#   include "initFluidVolume.H"
#   include "initLorentzForce.H"
#   include "initContinuityErrs.H"

    using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << nl << "Starting time loop" << nl << endl;

    while (runTime.run())
    {
        // Make the fluxes relative
        fvc::makeRelative(phi,U);

#       include "readGravitationalAcceleration.H"
#       include "readControls.H"

#       include "updateCourantNumber.H"
#       include "updateDeltaT.H"

        // Make the fluxes absolute
        fvc::makeAbsolute(phi,U);

        runTime++;

        Info << "Time = " << runTime.timeName() << nl << endl;

        interface.moveMeshPointsForOldTrackedSurfDisplacement();

        interface.updateDisplacementDirections();

        interface.predictPoints();

        // --- PIMPLE loop
        for (int outerCorr=0; outerCorr<max(nOuterCorr,nMinOuterCorr); outerCorr++)
        {
            // Initialize block matrix and set ldu to "zero"
#           include "initBlockMatrix.H"

            // Correct fixedGradient boundary conditions
#           include "updateLorentzForcePressureCorr.H"

            // Update interface bc
            interface.updateBoundaryConditions();

            // Make the fluxes relative
            fvc::makeRelative(phi,U);

#           include "updateCourantNumber.H"

#           include "UBlockEqn.H"

#           include "pBlockEqn.H"

#           include "UpBlockEqn.H"

            interface.correctPoints();

#           include "updateSwirl.H"
#           include "updateTurbulence.H"
#           include "updateVolErrorCheck.H"
        }

#       include "updateInterfaceForceCheck.H"
#       include "updateLorentzForce.H"
#       include "updateFluidVolume.H"

//         if (runTime.outputTime()) interface.writeVolA();

        runTime.write();

        Info << "ExecutionTime = "
             << scalar(runTime.elapsedCpuTime())
             << " s" << nl << endl;
    }

    Info << "End" << nl << endl;

    delete pBfacPtr; pBfacPtr = NULL;

    return(0);
}

// ************************************************************************* //
