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
    interTrackManager

Description
    ...

SourceFiles
    interTrackManager.C
    interTrackManagerDefaultRegionStorage.C

\*---------------------------------------------------------------------------*/

#ifndef interTrackManager_H
#define interTrackManager_H

#include "solverManager.H"

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "interTrackControl.H"
#include "twoPhaseMixture.H"
#include "turbulenceModel.H"
#include "trackedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class interTrackManager Declaration
\*---------------------------------------------------------------------------*/

class interTrackManager
:
    public solverManager<dynamicFvMesh>
{
public:

    // Public Nested Classes

// TODO: Move DefaultRegion declaration into Regions?
        //- Default region
        smRegion_new(interTrackManager, DefaultRegion);

// TODO: Probably its better to use a hash table for Settings?
// TODO: Read this settings based on Region class? -> init(), read(), ...
            //- Settings
            public: struct Settings
            {
                Switch UpCoupled;
                Switch UpPredictorForce;
            };
            smRegion_addSettings();

            //- Control
            smRegion_addControl_mesh(interTrackControl);

            //- Storage
            smStorage_new(interTrackManager, Storage);

                //- Gravitational acceleration
                smStorage_add(uniformDimensionedVectorField, g);

                //- Pressure
                smStorage_add(volScalarField, p);

                //- Velocity
                smStorage_add(volVectorField, U);

                //- Face flux
                smStorage_add(surfaceScalarField, phi);

                //- Volume force
                smStorage_add(volVectorField, F);

                //- Density
                smStorage_add(volScalarField, rho);

                //- Dynamic viscosity
                smStorage_add(volScalarField, mu);

                //- Fluid indicator
                smStorage_add(volScalarField, fluidIndicator);

                //- Transport model
                smStorage_add(twoPhaseMixture, transport);

                //- Turbulence model
                smStorage_add(incompressible::turbulenceModel, turbulence);

                //- Tracked surface
                smStorage_add(trackedSurface, interface);

            smStorage_end();
            smRegion_addStorage_manager();

        smRegion_end();

        //- Regions
        smRegions_new(interTrackManager, Regions);

            //- Default region
            smRegions_addDefaultRegion_manager(DefaultRegion);

        smRegions_end();

        //- Virtual regions function
        smBaseNested_addVirtualItem_P1(Regions, regions, *this);


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        interTrackManager(const interTrackManager&);

        //- Disallow default bitwise assignment
        void operator=(const interTrackManager&);

        // Evolution

            //- Set Courant Number
            virtual bool setCoNum(scalar& CoNum) const;

            //- Init sequence
            virtual void init() const;

            //- Read sequence
            virtual void read() const;

            //- Finalize sequence
            virtual void finalize() const;

            //- Pre write sequence
            virtual void writePre() const;

            //- Post write sequence
            virtual void writePost() const;

public:

    // Runtime type information
    TypeName("interTrackManager");

    // Constructor

        interTrackManager
        (
            const argList& args,
            Time& time,
            dynamicFvMesh& mesh,
            const word& prefix = "interTrack",
            const bool& master = true
        );

    // Destructor

        virtual ~interTrackManager() {};

    // Member Functions
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "interTrackManager.C"
#   include "interTrackManagerDefaultRegionStorage.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
