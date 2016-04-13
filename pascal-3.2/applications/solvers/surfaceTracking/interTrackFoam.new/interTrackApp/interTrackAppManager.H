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
    Foam::interTrackApp::Manager

Description
    ...

SourceFiles
    interTrackAppManager.C
    interTrackAppManagerRegion_DEFAULT.C

\*---------------------------------------------------------------------------*/

#ifndef interTrackAppManager_H
#define interTrackAppManager_H

#include "interTrackAppInclude.H"
#include "interTrackAppRegion.H"
#include "interTrackAppControl.H"

#include "solverManager.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace interTrackApp
{

/*---------------------------------------------------------------------------*\
                   Class interTrackApp::Manager Declaration
\*---------------------------------------------------------------------------*/

class Manager
:
    public solverManager<dynamicFvMesh, Region::SIZE>
{
public:

    // Public typedefs and forward declarations

        typedef interTrackApp::Control Control;
        class Settings;
        class Storage;

    // Solver Manager Construction

        SM_GLOBALCONTROL_ADD_P1(Control, mesh())


        SM_GLOBALSETTINGS_NEW();

            //- Debug switch
            SM_SETTINGS_MEMBER_ADD(Switch, debug);

        SM_GLOBALSETTINGS_ADD();


        SM_GLOBALSTORAGE_NEW();
        SM_GLOBALSTORAGE_ADD();


        SM_REGION_NEW(DEFAULT);

            SM_REGIONCONTROL_ADD_GLOBALCONTROL();

            SM_REGIONSETTINGS_NEW();

                //- Volume force switch
                SM_SETTINGS_MEMBER_ADD(Switch, volumeForce);

                //- Turbulence switch
                SM_SETTINGS_MEMBER_ADD(Switch, turbulence);

                //- Coupled solution of U and p
                SM_SETTINGS_MEMBER_ADD(Switch, UpCoupled);

                //- Directly inject volume force F in U-equation as
                //  source term
                //
                //  * enabled: UEqn -= F
                //
                //  * disabled: PhiHbyA += (rAU*F)_f
                //
                SM_SETTINGS_MEMBER_ADD(Switch, UEqnVolumeForce);

                //- Use fixed flux-pressure to correct pressure p at
                //  fixedGradient patches
                // TODO: Do not use fixedGradient
                //
                //  * enabled: dp/dn = (PhiHbyA -(U&Sf))/(magSf*rAUf)
                //
                //  * disabled: dp/dn = F*n
                //
                SM_SETTINGS_MEMBER_ADD(Switch, snGradpFromFlux);

            SM_REGIONSETTINGS_ADD();

            SM_REGIONSTORAGE_NEW();

                //- Gravitational acceleration
                SM_REGIONITEM_ADD(uniformDimensionedVectorField, g);

                //- Pressure
                SM_REGIONITEM_ADD(volScalarField, p);

                //- Velocity
                SM_REGIONITEM_ADD(volVectorField, U);

                //- Velocity/Pressure
                SM_REGIONITEM_ADD(volVector4Field, Up);

                //- Face flux
                SM_REGIONITEM_ADD(surfaceScalarField, phi);

                //- Volume force
                SM_REGIONITEM_ADD(volVectorField, F);

                //- Density
                SM_REGIONITEM_ADD(volScalarField, rho);

                //- Dynamic viscosity
                SM_REGIONITEM_ADD(volScalarField, mu);

                //- Fluid indicator
                SM_REGIONITEM_ADD(volScalarField, fluidIndicator);

                //- Transport model
                SM_REGIONITEM_ADD(twoPhaseMixture, transport);

                //- Turbulence model
                SM_REGIONITEM_ADD(incompressible::turbulenceModel, turbulence);

                //- Tracked surface
                SM_REGIONITEM_ADD(trackedSurface, interface);

            SM_REGIONSTORAGE_ADD();

        SM_REGION_END();


        SM_REGIONS_NEW();

            SM_REGIONS_REGION_ADD(DEFAULT, mesh());

        SM_REGIONS_ADD();


private:

    // Private Data

        // Regions

            //- Region mesh names
            wordList regionNames_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        Manager(const Manager&);

        //- Disallow default bitwise assignment
        void operator=(const Manager&);

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
    TypeName("interTrackApp::Manager");

    // Constructor

        Manager
        (
            const argList& args,
            Time& time,
            dynamicFvMesh& mesh,
            const word& prefix = "interTrack",
            const bool& master = true
        );

    // Destructor

        virtual ~Manager();

    // Member Functions

        // Regions

            //- Region mesh names
            inline const wordList& regionNames() const
            {
                return regionNames_;
            };

            //- Region mesh name for given region label
            inline word regionName(const label& regionI) const
            {
                return regionNames_[regionI];
            };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interTrackApp

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "interTrackAppManager.C"
#   include "interTrackAppManagerRegion_DEFAULT.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
