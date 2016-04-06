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
    eddyCurrentManager

Description
    ...

SourceFiles
    eddyCurrentManager.C
    eddyCurrentManagerDefaultRegion.C
    eddyCurrentManagerConductorRegion.C

\*---------------------------------------------------------------------------*/

#ifndef eddyCurrentManager_H
#define eddyCurrentManager_H

#include "solverManager.H"

#include "eddyCurrentControl.H"

#include "fvCFD.H"
#include "regionModelling.H"
#include "fvBlockMatrix.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class eddyCurrentManager Declaration
\*---------------------------------------------------------------------------*/

class eddyCurrentManager
:
    public solverManagerCore<regionFvMesh>
{
public:

    // Public data types

        SM_GLOBALCONTROL_ADD_P1(eddyCurrentControl, mesh()[polyMesh::defaultRegion])

        SM_GLOBALSETTINGS_NEW(eddyCurrentManager);

            //- Debug switch
            SM_SETTINGS_MEMBER_ADD(Switch, debug);

        SM_GLOBALSETTINGS_ADD();


        SM_GLOBALSTORAGE_NEW(eddyCurrentManager);

            //- Harmonic field frequency
            SM_STORAGE_MEMBER_ADD(uniformDimensionedScalarField, f0);

            //- Angular harmonic field frequency
            SM_STORAGE_MEMBER_ADD(uniformDimensionedScalarField, omega0);

            //- Electrical conductivity
            SM_STORAGE_MEMBER_ADD(regionVolScalarField, sigma);

            //- Magnetic vector potential
            SM_STORAGE_MEMBER_ADD(regionVolVectorField, ARe);
            SM_STORAGE_MEMBER_ADD(regionVolVectorField, AIm);

            //- Gradient of electric scalar potential
            SM_STORAGE_MEMBER_ADD(regionVolVectorField, VReGrad);
            SM_STORAGE_MEMBER_ADD(regionVolVectorField, VImGrad);

            //- Lorentz Force (time-averaged)
            SM_STORAGE_MEMBER_ADD(regionVolVectorField, FL);

            //- Magnetic pressure (time-averaged)
            SM_STORAGE_MEMBER_ADD(regionVolScalarField, pB);

        SM_GLOBALSTORAGE_ADD();


        SM_DEFAULTREGION_NEW(eddyCurrentManager);

            SM_SETTINGS_NEW(eddyCurrentManager);
            SM_SETTINGS_ADD();

            SM_STORAGE_NEW(eddyCurrentManager);

                // Source current density
                SM_STORAGE_MEMBER_ADD(volVectorField, j0Re);
                SM_STORAGE_MEMBER_ADD(volVectorField, j0Im);

                // Current density
                SM_STORAGE_MEMBER_ADD(volVectorField, jRe);
                SM_STORAGE_MEMBER_ADD(volVectorField, jIm);

                // Magnetic field density
                SM_STORAGE_MEMBER_ADD(volVectorField, BRe);
                SM_STORAGE_MEMBER_ADD(volVectorField, BIm);

            SM_STORAGE_ADD();

        SM_REGION_END();


        SM_REGION_NEW(eddyCurrentManager, ConductorRegion);

            SM_SETTINGS_NEW(eddyCurrentManager);
            SM_SETTINGS_ADD();

            SM_STORAGE_NEW(eddyCurrentManager);

                //- Electric scalar potential
                SM_STORAGE_MEMBER_ADD(volScalarField, VRe);
                SM_STORAGE_MEMBER_ADD(volScalarField, VIm);

                //- Electric scalar potential average
                SM_STORAGE_MEMBER_ADD(uniformDimensionedScalarField, VReAverage);
                SM_STORAGE_MEMBER_ADD(uniformDimensionedScalarField, VImAverage);

            SM_STORAGE_ADD();

        SM_REGION_END();


        SM_REGIONS_NEW(eddyCurrentManager);

            SM_REGIONS_DEFAULTREGION_ADD(mesh()[polyMesh::defaultRegion]);

// TODO: Read name during construction
            SM_REGIONS_REGION_ADD(ConductorRegion, conductorRegion, mesh()["region_conductor"]);

        SM_REGIONS_ADD();


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        eddyCurrentManager(const eddyCurrentManager&);

        //- Disallow default bitwise assignment
        void operator=(const eddyCurrentManager&);

        // Evolution

            //- Next sequence
            virtual void next() const;

            //- Write sequence
            virtual void write() const;

            //- Finalize sequence
            virtual void finalize() const;


public:

    // Runtime type information
    TypeName("eddyCurrentManager");

    // Constructor

        eddyCurrentManager
        (
            const argList& args,
            Time& time,
            regionFvMesh& mesh,
            const word& prefix = "eddyCurrent",
            const bool& master = true
        );

    // Destructor

        virtual ~eddyCurrentManager() {};

    // Member Functions
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "eddyCurrentManager.C"
#   include "eddyCurrentManagerDefaultRegion.C"
#   include "eddyCurrentManagerConductorRegion.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
