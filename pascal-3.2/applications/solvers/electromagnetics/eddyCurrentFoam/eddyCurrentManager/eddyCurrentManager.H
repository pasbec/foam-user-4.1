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
    eddyCurrentManagerRegionDEFAULT.C
    eddyCurrentManagerRegionCONDUCTOR.C

\*---------------------------------------------------------------------------*/

#ifndef eddyCurrentManager_H
#define eddyCurrentManager_H

#include "solverManager.H"

#include "eddyCurrentControl.H"

#include "fvCFD.H"
#include "regionModelling.H"
#include "fvBlockMatrix.H"
#include "physicalConstants.H"
#include "uniformDimensionedFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace solverManagerRegion
{
    enum Labels
    {
        DEFAULT,
        CONDUCTOR
    };

    static const label SIZE = 2;
}

/*---------------------------------------------------------------------------*\
                       Class eddyCurrentManager Declaration
\*---------------------------------------------------------------------------*/

class eddyCurrentManager
:
    public solverManager<regionFvMesh, solverManagerRegion::SIZE>
{
public:

    // Solver Manager Construction

        SM_HEADER_ADD(eddyCurrentManager)


        SM_GLOBALCONTROL_ADD_P1(eddyCurrentControl, mesh()[solverManagerRegion::DEFAULT])


        SM_GLOBALSETTINGS_NEW();

            //- Debug switch
            SM_SETTINGS_MEMBER_ADD(Switch, debug);

        SM_GLOBALSETTINGS_ADD();


        SM_GLOBALSTORAGE_NEW();

            //- Harmonic field frequency
            SM_GLOBALITEM_ADD(uniformDimensionedScalarField, f0);

            //- Angular harmonic field frequency
            SM_GLOBALITEM_ADD(uniformDimensionedScalarField, omega0);

            //- Electrical conductivity
            SM_GLOBALITEM_ADD(regionVolScalarField, sigma);

            //- Magnetic vector potential
            SM_GLOBALITEM_ADD(regionVolVectorField, ARe);
            SM_GLOBALITEM_ADD(regionVolVectorField, AIm);

            //- Gradient of electric scalar potential
            SM_GLOBALITEM_ADD(regionVolVectorField, VReGrad);
            SM_GLOBALITEM_ADD(regionVolVectorField, VImGrad);

            //- Lorentz Force (time-averaged)
            SM_GLOBALITEM_ADD(regionVolVectorField, FL);

            //- Magnetic pressure (time-averaged)
            SM_GLOBALITEM_ADD(regionVolScalarField, pB);

        SM_GLOBALSTORAGE_ADD();


        SM_REGION_NEW(DEFAULT);

            SM_REGIONCONTROL_ADD_GLOBALCONTROL();

            SM_REGIONSETTINGS_NEW();
            SM_REGIONSETTINGS_ADD();

            SM_REGIONSTORAGE_NEW();

                // Source current density
                SM_REGIONITEM_ADD(volVectorField, j0Re);
                SM_REGIONITEM_ADD(volVectorField, j0Im);

                // Current density
                SM_REGIONITEM_ADD(volVectorField, jRe);
                SM_REGIONITEM_ADD(volVectorField, jIm);

                // Magnetic field density
                SM_REGIONITEM_ADD(volVectorField, BRe);
                SM_REGIONITEM_ADD(volVectorField, BIm);

            SM_REGIONSTORAGE_ADD();

        SM_REGION_END();


        SM_REGION_NEW(CONDUCTOR);

            SM_REGIONCONTROL_ADD_GLOBALCONTROL();

            SM_REGIONSETTINGS_NEW();
            SM_REGIONSETTINGS_ADD();

            SM_REGIONSTORAGE_NEW();

                //- Electric scalar potential
                SM_REGIONITEM_ADD(volScalarField, VRe);
                SM_REGIONITEM_ADD(volScalarField, VIm);

                //- Electric scalar potential average
                SM_REGIONITEM_ADD(uniformDimensionedScalarField, VReAverage);
                SM_REGIONITEM_ADD(uniformDimensionedScalarField, VImAverage);

            SM_REGIONSTORAGE_ADD();

        SM_REGION_END();


        SM_REGIONS_NEW();

            SM_REGIONS_REGION_ADD(DEFAULT, mesh()[solverManagerRegion::DEFAULT]);

            SM_REGIONS_REGION_ADD(CONDUCTOR, mesh()[solverManagerRegion::CONDUCTOR]);

        SM_REGIONS_ADD();


private:

    // Private Data

        // Regions

            //- Region mesh names
            wordList regionNames_;

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

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "eddyCurrentManager.C"
#   include "eddyCurrentManagerRegionDEFAULT.C"
#   include "eddyCurrentManagerRegionCONDUCTOR.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
