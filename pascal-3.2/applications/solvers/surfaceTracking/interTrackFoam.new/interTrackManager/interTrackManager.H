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

        //- Default region
        solverManagerRegion_new(interTrackManager, DefaultRegion);

            //- Control
            solverManagerRegion_addControl(interTrackControl, this->mesh());

            //- Storage class
            solverManagerStorage_new(interTrackManager, Storage);

                //- Gravitational acceleration
                solverManagerStorage_add(uniformDimensionedVectorField, g);

                //- Pressure
                solverManagerStorage_add(volScalarField, p);

                //- Velocity
                solverManagerStorage_add(volVectorField, U);

                //- Face flux
                solverManagerStorage_add(surfaceScalarField, phi);

                //- Volume force
                solverManagerStorage_add(volVectorField, F);

                //- Density
                solverManagerStorage_add(volScalarField, rho);

                //- Dynamic viscosity
                solverManagerStorage_add(volScalarField, mu);

                //- Fluid indicator
                solverManagerStorage_add(volScalarField, fluidIndicator);

                //- Transport model
                solverManagerStorage_add(twoPhaseMixture, transport);

                //- Turbulence model
                solverManagerStorage_add(incompressible::turbulenceModel, turbulence);

                //- Tracked surface
                solverManagerStorage_add(trackedSurface, interface);

            solverManagerStorage_end();

            solverManagerRegion_addStorage();

        solverManagerRegion_end();

        //- Regions
        solverManagerRegions_new(interTrackManager, Regions);

            //- Default region
            solverManagerRegions_addDefaultRegion(DefaultRegion);

        solverManagerRegions_end();

    // Add virtual regions function
    solverManager_addVirtual(Regions, regions);


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        interTrackManager(const interTrackManager&);

        //- Disallow default bitwise assignment
        void operator=(const interTrackManager&);

        // Evolution

            //- Set Courant Number
            virtual bool setCoNum(scalar& CoNum) const;

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

        //- Read
        virtual void init() const;
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
