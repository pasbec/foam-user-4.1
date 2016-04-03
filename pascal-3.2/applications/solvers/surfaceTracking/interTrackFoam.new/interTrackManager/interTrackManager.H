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
    interTrackManagerDefaultRegion.C

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

        //- Settings
        smSettings_new(interTrackManager, Settings);

            mutable Switch debug;

        smSettings_end();
        smManager_addSettings_this();


        //- Default region
        smRegion_new(interTrackManager, DefaultRegion);

            //- Control
            smRegion_addControl_mesh(interTrackControl);

            //- Settings
            smSettings_new(interTrackManager, Settings);

                mutable Switch UpCoupled;
                mutable Switch UpDirectForce;

            smSettings_end();
            smRegion_addSettings_manager();

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
        smManager_addRegions_this();


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        interTrackManager(const interTrackManager&);

        //- Disallow default bitwise assignment
        void operator=(const interTrackManager&);

        // Evolution

            //- Set Courant Number
            virtual bool setCoNum(scalar& CoNum) const;

            //- Next sequence
            virtual void next() const;

            //- Write sequence
            virtual void write() const;

            //- Finalize sequence
            virtual void finalize() const;


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
#   include "interTrackManagerDefaultRegion.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
