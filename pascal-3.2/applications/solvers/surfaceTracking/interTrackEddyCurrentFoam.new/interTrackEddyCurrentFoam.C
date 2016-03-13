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
#include "blockMatrixTools.H" // TODO
#include "interTrackEddyCurrentControl.H"

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
#   include "addTrackedSurfacePrefixOption.H"

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createRegionDynamicFvMesh.H"

#   include "interTrackEddyCurrentRegionFields.H"
#   include "eddyCurrentConductorFields.H"

    // TODO: From here!!!

    // TODO: Solver debug flag
    bool debug = true;

#   include "initRegionCourantNo.H"

#   include "createBaseFields.H"
#   include "readBaseControls.H"
#   include "initBaseMagnetic.H"

#   include "createFluidFields.H"
#   include "readFluidControls.H"
#   include "initFluidContinuityErrs.H"
#   include "initFluidVolume.H"

#ifndef namespaceFoam
#define namespaceFoam
    using namespace Foam;
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    control.msg().startTimeLoop();

    while (runTime.run())
    {
        control.msg().timeStep();

// TODO [High]: Do CourantNo reset in regionControl
        // Reset global Courant No
#       include "initRegionCourantNo.H"

        // ==================================================================//
        // Re-read variable controls and adjust time step
        // ==================================================================//

#       include "readBaseControls.H"
        {
            setRegionScope(control.base());

#           include "updateRegionCourantNo.H"
        }

#       include "readFluidControls.H"
        {
            setRegionScope(control.fluidA());

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

        control.msg().time();

// // TODO TEST START
//         const surfaceScalarField& faceDiffusivityFluid =
//             mesh_[control.fluidA()].lookupObject<surfaceScalarField>("faceDiffusivity");
//
//         volScalarField diffusivityFluid =
//             volScalarField
//             (
//                 IOobject
//                 (
//                     "diffusivity",
//                     runTime.timeName(),
//                     mesh_[control.fluidA()],
//                     IOobject::NO_READ,
//                     IOobject::AUTO_WRITE
//                 ),
//                 fvc::average(faceDiffusivityFluid)()
//             );
//
//         const surfaceScalarField& faceDiffusivityDynamic =
//             mesh_[control.dynamic()].lookupObject<surfaceScalarField>("faceDiffusivity");
//
//         volScalarField diffusivityDynamic =
//             volScalarField
//             (
//                 IOobject
//                 (
//                     "diffusivity",
//                     runTime.timeName(),
//                     mesh_[control.dynamic()],
//                     IOobject::NO_READ,
//                     IOobject::AUTO_WRITE
//                 ),
//                 fvc::average(faceDiffusivityDynamic)()
//             );
// // TODO TEST END

        // ==================================================================//
        // Move mesh of fluid region and do prediction step
        // ==================================================================//

        {
            setRegionScope(control.fluidA());

            // Move and update mesh of fluid region
            interface.moveMeshPointsForOldTrackedSurfDisplacement();
            interface.updateDisplacementDirections();
            interface.updateControlPointsPosition();

// // TODO TEST: Make sure phi is correct
//             // Make sure phi is correct at interface
//             phi.boundaryField()[interface.aPatchID()] =
//                 fvc::meshPhi(U)().boundaryField()[interface.aPatchID()];

//             // Prediction step for interface points
//             interface.predictPoints();
        }

//         // TODO
//         runTime.writeNow();

        // ==================================================================//
        // Move mesh of dynamic region
        // ==================================================================//

// TODO: WARNING - THEORETICALLY,
//                 We have to apply point displacement during predictPoints for
//                 the dynamic region here. But predictPoints is USELESS!!!
//                 Currently it does NOTHING as (phi == meshPhi)@trackedSurface


        // ==================================================================//
        // Check whether magnetic update is due
        // ==================================================================//

        {
            // If magnetic update is due, mfUpdate will be true
#           include "updateBaseMagneticCheck.H"
        }

        // ==================================================================//
        // Update mesh points of base region
        // ==================================================================//

        if (mfUpdate)
        {
            // Create new point field for base region
            // with current point positions of fluid region
            pointField newPoints = mesh_.rmap(control.fluidA());

            // Replace point positions of dynamic region in
            // new point field for base region
            mesh_.rmap(newPoints,control.dynamic());

            // Update base mesh
            mesh_[control.base()].movePoints(newPoints);
            mesh_[control.base()].update();
        }

        // ==================================================================//
        // Update mesh points of conductor region if magnetic update is due
        // ==================================================================//

        if (mfUpdate)
        {
            // Create new point field for conductor region
            // with current points of base region
            pointField newPoints = mesh_.map(control.conductor());

            // Move and update mesh of conductor region
            mesh_[control.conductor()].movePoints(newPoints);
            mesh_[control.conductor()].update();
        }

        // ==================================================================//
        // Solve for magnetic and derived fields in base/conductor region
        // ==================================================================//

        if (mfUpdate)
        {
            setControlScope(eddyCurrent);

#           include "eddyCurrentAVloop.H"

#           include "eddyCurrentDerivedFields.H"
        }

        // ==================================================================//
        // Map Lorentz-force and magnetic pressure fields from base region
        // to fluid region and extrapolate their boundary values
        // ==================================================================//

        if (mfUpdate)
        {
//             FL_.mapExtrapolate(control.fluidA());
//             pB_.mapExtrapolate(control.fluidA());
// TODO TEST FIXME Extrapolation sometimes fails
            FL_.mapCopyInternal(control.fluidA());
            pB_.mapCopyInternal(control.fluidA());
        }

        // ==================================================================//
        // Process magnetic-field related fields in fluid region
        // and correct gradient of fluid pressure at corresponding boundaries
        // ==================================================================//

        if (mfUpdate)
        {
            setRegionScope(control.fluidA());

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
            setRegionScope(control.fluidA());

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
                control.fluidA(),
                control.dynamic()
            );

            // Grab current points of dynamic region as new point field
            pointField newPoints = mesh_[control.dynamic()].points();

// TODO/FIXME: Really necessary?
            // Correct points for 2D-motion of dynamic region
            twoDPointCorrector dynamicTwoDPointCorr(mesh_[control.dynamic()]);
            dynamicTwoDPointCorr.correctPoints(newPoints);

            // Use motionSolver to move mesh of dynamic region
            mesh_[control.dynamic()].movePoints(newPoints);
            mesh_[control.dynamic()].update();
        }

        // ==================================================================//
        // Finish time step and write
        // ==================================================================//

        runTime.write();

        if (debug)
        {
            if (runTime.outputTime())
            {
                interface.writeVTK();
                interface.writeVTKControlPoints();
                interface.writeVolA();
            }
        }

        control.msg().executionTime();
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    control.msg().end();

    return(0);
}

// ************************************************************************* //
