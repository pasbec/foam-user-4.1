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

        //- Global Settings
        SM_GLOBALSETTINGS_NEW(interTrackManager);

            //- Debug switch
            SM_GLOBALSETTINGS_MEMBER_ADD(Switch, debug);

        SM_GLOBALSETTINGS_ADD();


        //- Default region
        SM_DEFAULTREGION_NEW(interTrackManager);

            //- Control
            SM_REGION_CONTROL_ADD(interTrackControl);

            //- Settings
            SM_SETTINGS_NEW(interTrackManager);

                //- Coupled solution of U and p
                SM_SETTINGS_MEMBER_ADD(Switch, UpCoupled);

                //- Directly inject volume force F in U-equation
                SM_SETTINGS_MEMBER_ADD(Switch, UpDirectForce);

            SM_SETTINGS_ADD();

            //- Storage
            SM_STORAGE_NEW(interTrackManager);

                //- Gravitational acceleration
                SM_STORAGE_MEMBER_ADD(uniformDimensionedVectorField, g);

                //- Pressure
                SM_STORAGE_MEMBER_ADD(volScalarField, p);

                //- Velocity
                SM_STORAGE_MEMBER_ADD(volVectorField, U);

                //- Face flux
                SM_STORAGE_MEMBER_ADD(surfaceScalarField, phi);

                //- Volume force
                SM_STORAGE_MEMBER_ADD(volVectorField, F);

                //- Density
                SM_STORAGE_MEMBER_ADD(volScalarField, rho);

                //- Dynamic viscosity
                SM_STORAGE_MEMBER_ADD(volScalarField, mu);

                //- Fluid indicator
                SM_STORAGE_MEMBER_ADD(volScalarField, fluidIndicator);

                //- Transport model
                SM_STORAGE_MEMBER_ADD(twoPhaseMixture, transport);

                //- Turbulence model
                SM_STORAGE_MEMBER_ADD(incompressible::turbulenceModel, turbulence);

                //- Tracked surface
                SM_STORAGE_MEMBER_ADD(trackedSurface, interface);

            SM_STORAGE_ADD();

        SM_REGION_END();


        //- Regions
        SM_REGIONS_NEW(interTrackManager);

            //- Default region
            SM_REGIONS_DEFAULTREGION_ADD(this->mesh());

        SM_REGIONS_ADD();


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
