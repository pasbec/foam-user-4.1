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

Class
    trackedSurface

Description
    Implementation of surface (interface) tracking method
    based on moving mesh.

SourceFiles
    trackedSurface.C
    trackedSurfacePointDisplacement.C
    makeTrackedSurfaceData.C

\*---------------------------------------------------------------------------*/

#ifndef TrackedSurface_H
#define TrackedSurface_H


#include "fvCFD.H"
#include "IOdictionary.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "faCFD.H"
#include "IOpatchToPatchInterpolation.H"
#include "dynamicFvMesh.H"
#include "twoPhaseMixture.H"
#include "turbulenceModel.H"

#include "surfactantProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class trackedSurface Declaration
\*---------------------------------------------------------------------------*/

class trackedSurface
:
    public IOdictionary
{
    // Private data

        //- Reference to fvMesh
        dynamicFvMesh& mesh_;

        //- Reference to density field
        const volScalarField& rho_;

        //- Reference to velocity field
        volVectorField& U_;

        //- Reference to pressure field
        volScalarField& p_;

        // Reference to optional external pressure field (ine fluid only)
        const volScalarField& pExt_;

        //- Reference to fluid flux field
        const surfaceScalarField& phi_;

        //- Reference to gravitational acceleration
        const dimensionedVector g_;

        //- Reference to the transportModel
        const twoPhaseMixture& transport_;

        // Reference to auoPtr of turbulence model
        const autoPtr<incompressible::turbulenceModel>& turbulence_;

        //- Current time index for surface
        label curTimeIndex_;

        //- Are there exist two fluids in the surface model
        Switch twoFluids_;

        //- Is it surface points displacement direction
        //  parallel with surface point normals
        Switch normalMotionDir_;

        //- Free-surface points displacement direction
        //  if not normal motion direction
        vector motionDir_;

        //- Is it surface clean of surfactant
        Switch cleanInterface_;

        //- ID of the fvMesh patch which represents
        //  A side of the surface
        label aPatchID_;

        //- ID of the fvMesh patch which represents
        //  B side of the surface
        label bPatchID_;

        //- Density of fluid A
        dimensionedScalar rhoFluidA_;

        //- Density of fluid B
        dimensionedScalar rhoFluidB_;

        //- Viscosity of fluid A
        dimensionedScalar muFluidA_;

        //- Viscosity of fluid B
        dimensionedScalar muFluidB_;

        //- Surface tension for the clean surface
        dimensionedScalar cleanInterfaceSurfTension_;

        //- surface patches which do not move
        wordList fixedTrackedSurfacePatches_;

        //- surface patches for witch point normals
        //  must be corrected
        wordList pointNormalsCorrectionPatches_;

         //- Number of surface correctors
        label nFreeSurfCorr_;

        //- Interface smoothing at the begining of time step
        Switch smoothing_;

        //- Are we using clean tangential surface tension correction
        Switch noCleanTangentialSurfaceTensionCorrection_;

        //- Are we using pointDisplacementSM instead of least squares fit
        Switch motionByPointDisplacementSM_;

    // Demand-driven data

        //- Patch to patch interpolation object which deals with
        //  interpolation of flow variables between two surface
        //  (interface) patches
        mutable IOpatchToPatchInterpolation* interpolatorABPtr_;

        //- Patch to patch interpolation object which deals with
        //  interpolation of flow variables between two surface
        //  (interface) patches
        mutable IOpatchToPatchInterpolation* interpolatorBAPtr_;

        //- Points which are attached to the surface A side faces
        //  and which defines the surface shape
        mutable vectorIOField* controlPointsPtr_;

        //- Field which additionally determines
        //  the motion of surface points
        mutable labelList* motionPointsMaskPtr_;

        //- Displacement direction of surface points
        mutable vectorField* pointsDisplacementDirPtr_;

        //- Displacement direction of surface control points
        mutable vectorField* facesDisplacementDirPtr_;

        //- Total displacement of surface points
        //  in one time step
        mutable vectorIOField* totalDisplacementPtr_;

        //- Finite area mesh which consists
        //  surface side A faces
        mutable faMesh* aMeshPtr_;

        //- Free-surface velocity field
        mutable areaVectorField* UsPtr_;

        //- Free-surface fluid flux
        mutable edgeScalarField* phisPtr_;

        //- Free-surface surfactant concetration
        mutable areaScalarField* surfactConcPtr_;

        //- Surface tension field
        mutable areaScalarField* surfaceTensionPtr_;

        //- Surfactant properties
        mutable surfactantProperties* surfactantPtr_;

        //- Fluid indicator
        mutable volScalarField* fluidIndicatorPtr_;

        //- Effiective viscosity value field of fluid A
        mutable scalarField* muEffFluidAvalPtr_;

        //- Effiective viscosity value field of fluid B
        mutable scalarField* muEffFluidBvalPtr_;

    // Private Member Functions

        // Make demand-driven data

            void makeInterpolators();
            void makeControlPoints();
            void makeMotionPointsMask();
            void makeDirections();
            void makeTotalDisplacement();
            void readTotalDisplacement();
            void makeFaMesh() const;
            void makeUs() const;
            void makePhis();
            void makePhi();
            void makeDdtPhi();
            void makeSurfactConc() const;
            void makeSurfaceTension() const;
            void makeSurfactant() const;
            void makeFluidIndicator() const;
            void makeMuEffFluidAval() const;
            void makeMuEffFluidBval() const;

        //- Clear all demand-driven data
        void clearOut();

        //- Disallow default bitwise copy construct
        trackedSurface(const trackedSurface&);

        //- Disallow default bitwise assignment
        void operator=(const trackedSurface&);

        //- Move control ponts for deltaH and calculate
        //  surface points displacement
        //  for given new control points position
        tmp<vectorField> pointDisplacement(const scalarField& deltaH);

        //-
        tmp<vectorField> pointDisplacementLSF(const scalarField& deltaH);

        // Calc least sqare plane point and normal
        tmp<vectorField> lsPlanePointAndNormal
        (
            const vectorField& points,
            const vector& origin,
            const vector& axis
        ) const;

        //-
        tmp<vectorField> pointDisplacementSM(const scalarField& deltaH);

public:

    // Declare name of the class and it's debug switch
    ClassName("trackedSurface");


    // Constructors

    trackedSurface
    (
        dynamicFvMesh& fvMesh,
        const volScalarField& volrho,
        volVectorField& volU,
        volScalarField& volp,
        const volScalarField& volpExt,
        const surfaceScalarField& sfphi,
        const uniformDimensionedVectorField unig,
        const twoPhaseMixture& transportModel,
        const autoPtr<incompressible::turbulenceModel>& pTurbulenceModel
    );


    // Destructor

    ~trackedSurface();


    // Member Functions


        //- Return reference to fvMesh
        dynamicFvMesh& mesh()
        {
            return mesh_;
        }

        //- Return reference to fvMesh
        const dynamicFvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return reference to DB
        const Time& DB() const
        {
            return mesh_.time();
        }

        //- Return reference to velocity field
        const volScalarField& rho() const
        {
            return rho_;
        }

        //- Return reference to velocity field
        volVectorField& U()
        {
            return U_;
        }

        //- Return reference to velocity field
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return reference to pressure field
        volScalarField& p()
        {
            return p_;
        }

        //- Return reference to pressure field
        const volScalarField& p() const
        {
            return p_;
        }

        //- Return reference to optional external pressure field
        const volScalarField& pExt() const
        {
            return pExt_;
        }

        //- Return refernece to gravity value
        const dimensionedVector& g() const
        {
            return g_;
        }

        //- Return reference to transport model
        const twoPhaseMixture& transport() const
        {
            return transport_;
        }

        //- Return reference to turbulence model
        const autoPtr<incompressible::turbulenceModel>& turbulence() const
        {
            return turbulence_;
        }

        //- Are there exist two fluid
        //  in the surface model
        bool twoFluids() const
        {
            return twoFluids_;
        }

        //- Motion direction swithc
        bool normalMotionDir() const
        {
            return normalMotionDir_;
        }

        //- Is it interface clean
        bool cleanInterface() const
        {
            return cleanInterface_;
        }

        //- ID of patch which represent A side
        //  of surface
        const label& aPatchID() const
        {
            return aPatchID_;
        }

        //- ID of patch which represent B side
        //  of surface
        const label& bPatchID() const
        {
            return bPatchID_;
        }

        //- Density of fluid A
        const dimensionedScalar& rhoFluidA() const
        {
            return rhoFluidA_;
        }

        //- Densiti of fluid B
        const dimensionedScalar& rhoFluidB() const
        {
            return rhoFluidB_;
        }

        //- Viscosity of fluid A
        const dimensionedScalar& muFluidA() const
        {
            return muFluidA_;
        }

        //- Viscosity of fluid B
        const dimensionedScalar& muFluidB() const
        {
            return muFluidB_;
        }

        //- Surface tension for clean interface
        const dimensionedScalar& cleanInterfaceSurfTension() const
        {
            return cleanInterfaceSurfTension_;
        }

        //- Are we using clean tangential surface tension correction
        bool noCleanTangentialSurfaceTensionCorrection() const
        {
            return noCleanTangentialSurfaceTensionCorrection_;
        }

        //- Are we using pointDisplacementSM instead of least squares fit
        bool motionByPointDisplacementSM() const
        {
            return motionByPointDisplacementSM_;
        }

        //- Return reference to interpolator from surface side A to
        //  surface side B
        const IOpatchToPatchInterpolation& interpolatorAB();

        //- Return reference to interpolator from surface side B to
        //  surface side A
        const IOpatchToPatchInterpolation& interpolatorBA();

        //- Return control points
        vectorField& controlPoints();

        //- Return reference to motion points mask field
        labelList& motionPointsMask();

        //- Return reference to point displacement direction field
        vectorField& pointsDisplacementDir();

        //- Return reference to control points displacement direction field
        vectorField& facesDisplacementDir();

        //- Return total points displacement
        vectorField& totalDisplacement();

        //- Return reference to finite area mesh
        faMesh& aMesh();

        //- Return reference to finite area mesh
        const faMesh& aMesh() const;

        //- Return surface velocity field
        areaVectorField& Us();

        //- Return surface velocity field
        const areaVectorField& Us() const;

        //- Return surface fluid flux field
        edgeScalarField& Phis();

        //- Return surface surfactant concentration field
        areaScalarField& surfactantConcentration();

        //- Return surface surfactant concentration field
        const areaScalarField& surfactantConcentration() const;

        //- Return surface tension field
        areaScalarField& surfaceTension();

        //- Return surface tension field
        const areaScalarField& surfaceTension() const;

        //- Return surface tension gradient
        tmp<areaVectorField> surfaceTensionGrad() const;

        //- Return surfactant properties
        const surfactantProperties& surfactant() const;

        //- Return reference to fluid indicator
        const volScalarField& fluidIndicator();

        //- Return effiective viscosity value field of fluid A
        scalarField& muEffFluidAval();

        //- Return effiective viscosity value field of fluid A
        const scalarField& muEffFluidAval() const;

        //- Return effiective viscosity value field of fluid B
        scalarField& muEffFluidBval();

        //- Return effiective viscosity value field of fluid B
        const scalarField& muEffFluidBval() const;

        //- Update control points end displacement directions
        void updateDisplacementDirections();

        //- Predict interface points position
        bool predictPoints
        (
            const scalar relax = 1.0
        );

        //- Correct interface points position
        bool correctPoints
        (
            const scalar relax = 1.0
        );

        //- Move only surface points
        bool movePoints
        (
            const scalarField& interfacePHi,
            const scalar relax = 1.0
        );

        //- Move mesh points for old time surface points displacement
        bool moveMeshPointsForOldFreeSurfDisplacement();

        //- Move mesh points
        bool moveMeshPoints
        (
            const scalar relax = 1.0
        );

        //- Update surface boundary conditions
        void updateBoundaryConditions();

        //- Correct effective kinematic viscosity on surface
        void updateMuEff();

        //- Correct surface pressure
        void updatePressure();

        //- Correct surface velocity
        void updateVelocity();

        //- Update surface flux
        void updateSurfaceFlux();

        //- Update surface surfactant concentration
        void updateSurfactantConcentration();

        //- Correct surface velocity boundary conditions
        void correctUsBoundaryConditions();

        //- Calculate total pressure force
        vector totalPressureForce() const;

        //- Calculate total viscous force
        vector totalViscousForce() const;

        //- Calculate total surface tension force
        vector totalSurfaceTensionForce() const;

        //- Calculate surface undulation indicator
        tmp<scalarField> undulationIndicator();

        //- Smooth surface
        void smooth();

        //- Chech surface faces flatness
        tmp<scalarField> checkFaceFlatness();

        //- Initialize control points position
        void initializeControlPointsPosition();

        //- Maximal surface tension based Courant number
        scalar maxCourantNumber();

        //- Global continuity error of tracked surface side A
        scalarField continuityErrorA() const;

        //- Global continuity error of tracked surface side B
        scalarField continuityErrorB() const;

        //
        template<class Type>
        void
        writeVolFromArea
        (
            const label aPatchID,
            const fvMesh& vMesh,
            const GeometricField<Type, faPatchField, areaMesh>& af
        );

        //- Write selected fields as volFields of tracked surface side A
        void writeVolA();

        //- Write VTK trackedSurface mesh
        void writeVTK() const;

        //- Write VTK trackedSurface control points
        void writeVTKControlPoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
