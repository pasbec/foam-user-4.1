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
#include "blockMatrixTools.H"

#include "regionDynamicFvMesh.H"
#include "regionVolFields.H"
#include "regionPointMesh.H"
#include "regionPointFields.H"
#include "regionControl.H"

#include "fixedGradientFvPatchFields.H"

#include "fvcExtrapolate.H"

#include "transportModel.H"
#include "twoPhaseMixture.H"
#include "turbulenceModel.H"

#include "trackedSurface.H"

#include "twoDPointCorrector.H"

// TODO [High]: Clearification:
//              - Why can we see velocity/pressure oscillations along
//                the interface at wall contact line at start?
//              - Why does Alberto-method result in different flow?
//              - Why is Alberto-method much more unstable?
//              - Why is "limited 0.5" for laplcian effecting
//                Alberto-method dramatically?

// TODO [High]: Structure/Implementation:
//              - Fix update of mesh points of base region (last main step).
//                The inner points of the base region get out of sync if
//                the fluid volume is corrected in the curent time step!
//              - Implement control class?
//              - Get rid of warning: unused variable for mesh, FL and pB
//              - Fix pressure relaxation for Alberto-method
//                (do we really need to store a field likeUcorrByrAU?)
//              - Info about main steps
//              - Find proper (better) solover name and rename it

// TODO [High]: Speed-Up/Memory:
//              - Move mesh only on mfUpdate?

// TODO [High]: Automatic remeshing:
//              - How to get an option/algorithm to redistribute all
//                mesh points in z-direction after some time?

// TODO [High]: Volume correction:
//              - Better way to realize this?
//              - Boundary points?
//                (If mass loss gets too large, boundary points might
//                mis-align during correction)

// TODO [Medium]: Add plausability checks:
//                 - interface.twoFluids()
//                 - patchField types and dimensions of all fields
//                   (especially calculated for FL and pB)
//                 - ...

// TODO [Low]: twoFluids:
//             - Make sure interface.twoFluids() must not be set?
//             - Fix all stuff with twoFluids
//             - Volume correction with twoFluids?

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Solver debug flag
    bool debug = true;

#   include "setRootCase.H"
#   include "createTime.H"

#   include "createRegionMesh.H"
#   include "createRegionFields.H"
#   include "createRegionControl.H"
#   include "initRegionCourantNo.H"

#   include "createBaseFields.H"
#   include "readBaseControls.H"
#   include "initBaseMagnetic.H"

#   include "createFluidFields.H"
#   include "readFluidControls.H"
#   include "initFluidContinuityErrs.H"
#   include "initFluidVolume.H"

    using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << nl << "Starting time loop" << nl << endl;

    while (runTime.run())
    {
        Info << "Time step = " << runTime.timeIndex() << nl << endl;

// TODO [High]: Do CourantNo reset in regionControl
        // Reset global Courant No
#       include "initRegionCourantNo.H"

        // ==================================================================//
        // Re-read variable controls and adjust time step
        // ==================================================================//

#       include "readBaseControls.H"
        {
#           include "setBaseFields.H"

#           include "updateRegionCourantNo.H"
        }

#       include "readFluidControls.H"
        {
#           include "setFluidFields.H"

#           include "updateFluidGravitationalAcceleration.H"
#           include "updateFluidCourantNumber.H"

#           include "updateRegionCourantNo.H"
        }

        {
// TODO [High]: Controls of dynamic region?
// TODO [High]: Mesh Courant Number of dynamic region?
        }

// TODO [High]: Global deltaT setting necessary?
// TODO [High]: Global deltaT setting in regionControl
#       include "updateRegionDeltaT.H"

        runTime++;

        Info << "Time = " << runTime.value() << nl << endl;

        // ==================================================================//
        // Move mesh of fluid region and do prediction step
        // ==================================================================//

        {
            // Move and update mesh of fluid region
            interface.moveMeshPointsForOldTrackedSurfDisplacement();
            interface.updateDisplacementDirections();

            // Prediction step for interface points
            interface.predictPoints();
        }

        // ==================================================================//
        // Check whether magnetic update is due
        // ==================================================================//

        {
            // If magnetic update is due, mfUpdate will be true
#           include "updateBaseMagneticCheck.H"
        }

        // ==================================================================//
        // Update mesh points of base region if magnetic update is due
        // ==================================================================//

        if (mfUpdate)
        {
            // Create new point field for base region
            // with current points of fluid region
            pointField newPoints = mesh_.rmap(fluidRegionID);

            // Move and update mesh of base region
            mesh_[defaultRegionID].movePoints(newPoints);
            mesh_[defaultRegionID].update();
        }

        // ==================================================================//
        // Solve for magnetic and derived fields in base region if due
        // ==================================================================//

        if (mfUpdate)
        {
#           include "setBaseFields.H"
#           include "solveBase.H"
        }

        // ==================================================================//
        // Map region fields from base region to fluid region
        // ==================================================================//

        if (mfUpdate)
        {
            FL_.map(fluidRegionID);
            pB_.map(fluidRegionID);
        }

        // ==================================================================//
        // Extrapolate boundary values of mapped fields in fluid region
        // and correct gradient of fluid pressure at corresponding boundaries
        // ==================================================================//

        if (mfUpdate)
        {
#           include "setFluidFields.H"

            fvc::extrapolate(FL);
            fvc::extrapolate(pB);

            pB *= lorentzForceRotationalFactor;

            F == FL + fvc::grad(pB);
            F *= lorentzForceVolumeFactor;

            if (debug)
            {
                volVectorField& debug_Fgrad = *debug_FgradPtr;

                debug_Fgrad == FL - F;
            }

#           include "updateFluidPressureBoundaryConditions.H"
        }

        // ==================================================================//
        // Solve fluid flow in fluid region and move interface
        // ==================================================================//

        {
#           include "setFluidFields.H"
#           include "solveFluid.H"
        }

        // ==================================================================//
        // Move mesh of dynamic region
        // ==================================================================//

        {
            // Calculate mesh velocity at dynamic/fluid-interface
            // in dynamic region from current boundary displacement
            // in fluid region
            mesh_.patchMapMeshVelocityDirectMapped
            (
                fluidRegionID,
                dynamicRegionID
            );

            // Grab current points of dynamic region as new point field
            pointField newPoints = mesh_[dynamicRegionID].points();

            // Correct points for 2D-motion of dynamic region
            twoDPointCorrector dynamicTwoDPointCorr(mesh_[dynamicRegionID]);
            dynamicTwoDPointCorr.correctPoints(newPoints);

            // Use motionSolver to move mesh of dynamic region
            mesh_[dynamicRegionID].movePoints(newPoints);
            mesh_[dynamicRegionID].update();
        }

        // ==================================================================//
        // Update mesh points of base region
        // ==================================================================//

        {
            // Create new point field for base region
            // with current points of patch trackedSurface in fluid region
// TODO [High] : If the fluid volume is being corrected in the current time step
//               we may not only use interface points, but all fluid points instead!
//               Do we really get loose much performance with this?
//             pointField newPoints = mesh_.rmap(fluidRegionID,"trackedSurface");
            pointField newPoints = mesh_.rmap(fluidRegionID);

            // Replace point positions of dynamic region in
            // new point field for base region
            mesh_.rmap(newPoints,dynamicRegionID);

            // Update base mesh
            mesh_[defaultRegionID].movePoints(newPoints);
            mesh_[defaultRegionID].update();

            Info << endl;
        }

        // ==================================================================//
        // Finish time step and write
        // ==================================================================//

        if (debug)
        {
            if (runTime.outputTime())
            {
                interface.writeVolA();
            }
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
