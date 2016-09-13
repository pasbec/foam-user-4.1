/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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
    Foam::eddyCurrentApp::Manager

Description
    ...

SourceFiles
    eddyCurrentAppManager.C
    eddyCurrentAppManagerRegion_DEFAULT.C
    eddyCurrentAppManagerRegion_CONDUCTOR.C

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#ifndef eddyCurrentAppManager_H
#define eddyCurrentAppManager_H

#include "eddyCurrentAppInclude.H"
#include "eddyCurrentAppControl.H"

#include "solverManagers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace eddyCurrentApp
{

/*---------------------------------------------------------------------------*\
                   Class eddyCurrentApp::Manager Declaration
\*---------------------------------------------------------------------------*/

class Manager
:
    public regionFvSolverManager
{
public:

    // Public typedefs and forward declarations

        typedef eddyCurrentApp::Control Control;
        class Settings;
        class Storage;

    // Solver Manager Construction

        SM_GLOBALCONTROL_ADD_P1(Control, mesh());


        SM_GLOBALSETTINGS_NEW();

            //- Switch to set wether tolerance scaling is used
            SM_SETTINGS_MEMBER_ADD(Switch, tolScale);

            //- Switch to set wether Biot-Savart is used
            SM_SETTINGS_MEMBER_ADD(Switch, biotSavart);

            //- Lorentz force switch
            SM_SETTINGS_MEMBER_ADD(Switch, lorentzForce);

            //- Magnetic pressure switch
            SM_SETTINGS_MEMBER_ADD(Switch, magneticPressure);

            //- Joul heat switch
            SM_SETTINGS_MEMBER_ADD(Switch, jouleHeat);

        SM_GLOBALSETTINGS_ADD();


        SM_GLOBALSTORAGE_NEW();

            //- Harmonic field frequency
            SM_GLOBALITEM_ADD(uniformDimensionedScalarField, f0);

            //- Angular harmonic field frequency
            SM_GLOBALITEM_ADD(uniformDimensionedScalarField, omega0);

            //- Electrical conductivity
            SM_GLOBALITEM_ADD(regionVolScalarField, sigma);

            //- Source vector potential
            SM_GLOBALITEM_ADD(regionVolVectorField, A0Re);
            SM_GLOBALITEM_ADD(regionVolVectorField, A0Im);

            //- Magnetic vector potential
            SM_GLOBALITEM_ADD(regionVolVectorField, ARe);
            SM_GLOBALITEM_ADD(regionVolVectorField, AIm);

            //- TODO
            SM_GLOBALITEM_ADD(uniformDimensionedVectorField, Anormf);
            SM_GLOBALITEM_ADD(uniformDimensionedVectorField, Ascale);

            //- Gradient of electric scalar potential
            SM_GLOBALITEM_ADD(regionVolVectorField, VReGrad);
            SM_GLOBALITEM_ADD(regionVolVectorField, VImGrad);

            // Magnetic field density
            SM_GLOBALITEM_ADD(regionVolVectorField, BRe);
            SM_GLOBALITEM_ADD(regionVolVectorField, BIm);

            //- Lorentz Force (time-averaged)
            SM_GLOBALITEM_ADD(regionVolVectorField, FL);

            //- Magnetic pressure (time-averaged)
            SM_GLOBALITEM_ADD(regionVolScalarField, pB);

            //- Joule heating (time-averaged)
            SM_GLOBALITEM_ADD(regionVolScalarField, Q);

        SM_GLOBALSTORAGE_ADD();


        SM_REGION_NEW(DEFAULT);

            SM_REGIONCONTROL_ADD_GLOBALCONTROL();

            SM_REGIONSETTINGS_NEW();
            SM_REGIONSETTINGS_ADD();

            SM_REGIONSTORAGE_NEW();

                // Source current density
                SM_REGIONITEM_ADD(volVectorField, j0Re);
                SM_REGIONITEM_ADD(volVectorField, j0Im);

            SM_REGIONSTORAGE_ADD();

        SM_REGION_END();


        SM_REGION_NEW(CONDUCTOR);

            SM_REGIONCONTROL_ADD_GLOBALCONTROL();

            SM_REGIONSETTINGS_NEW();
            SM_REGIONSETTINGS_ADD();

            SM_REGIONSTORAGE_NEW();

                //- Electric scalar potential
                SM_REGIONITEM_ADD(surfaceScalarField, ddtARePhi);
                SM_REGIONITEM_ADD(surfaceScalarField, ddtAImPhi);

                //- Electric scalar potential
                SM_REGIONITEM_ADD(volScalarField, VRe);
                SM_REGIONITEM_ADD(volScalarField, VIm);

                // Current density
                SM_REGIONITEM_ADD(volVectorField, jRe);
                SM_REGIONITEM_ADD(volVectorField, jIm);

            SM_REGIONSTORAGE_ADD();

        SM_REGION_END();


        SM_REGIONS_NEW();

            SM_REGIONS_REGION_ADD(DEFAULT, mesh()[Region::DEFAULT]);

            SM_REGIONS_REGION_ADD(CONDUCTOR, mesh()[Region::CONDUCTOR]);

        SM_REGIONS_ADD();


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        Manager(const Manager&);

        //- Disallow default bitwise assignment
        void operator=(const Manager&);


public:

    // Runtime type information
    TypeName("eddyCurrentApp::Manager");

    // Constructor

        Manager
        (
            const argList& args,
            Time& time,
            regionFvMesh& mesh,
            bool master = true,
            const word& prefix = "eddyCurrent"
        );

    // Destructor

        virtual ~Manager();

    // Member Functions

        // Evolution

            //- Next sequence
            virtual void next() const;

            //- Write sequence
            virtual void write() const;

            //- Finalize sequence
            virtual void finalize() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace eddyCurrentApp

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "eddyCurrentAppManager.C"
#   include "eddyCurrentAppManagerRegion_DEFAULT.C"
#   include "eddyCurrentAppManagerRegion_CONDUCTOR.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
