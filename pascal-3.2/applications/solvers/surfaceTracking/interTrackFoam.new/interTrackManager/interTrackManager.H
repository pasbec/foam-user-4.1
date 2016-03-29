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
    interTrackManagerStorage.C

\*---------------------------------------------------------------------------*/

#ifndef interTrackManager_H
#define interTrackManager_H

#include "solverManager.H"

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "trackedSurface.H"

#include "pimpleControl.H"

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

        //- Storage
        struct storage
        :
            solverManager<dynamicFvMesh>::storage
        {
        private:

            // Private Data

                //- Reference to solver manager
                interTrackManager& manager_;

            // Private Member Functions

                //- Disallow default bitwise copy construct
                storage(const storage&);

                //- Disallow default bitwise assignment
                void operator=(const storage&);

                //- Return reference to solver manager
                interTrackManager& manager()
                {
                    return manager_;
                }

                //- Initialize data sequence
                void init();

            // Macro-generated data entries with accessors

                //- Pressure
                solverManagerStorage_add(volScalarField, p);

                //- Velocity
                solverManagerStorage_add(volVectorField, U);

                //- Face flux
                solverManagerStorage_add(surfaceScalarField, phi);

                //- Density
                solverManagerStorage_add(volScalarField, rho);

                //- Dynamic viscosity
                solverManagerStorage_add(volScalarField, mu);

                //- Fluid indicator
                solverManagerStorage_add(volScalarField, fluidIndicator);

                //- Tracked surface
                solverManagerStorage_add(trackedSurface, interface);

        public:

            // Constructor

                //- Construct from manager
                storage(interTrackManager& manager);
        };

private:

    // Private Data

        //- Storage
        mutable storage* storagePtr_;

        //- Solution control
        mutable pimpleControl pimple_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        interTrackManager(const interTrackManager&);

        //- Disallow default bitwise assignment
        void operator=(const interTrackManager&);

        // Evolution

            //- Time step setting routine
            virtual bool calcCoNum(scalar& CoNum) const;

            //- Post write
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
            const word& name = "interTrack",
            const bool& master = true
        );

    // Destructor

        virtual ~interTrackManager()
        {
            if (this->master()) delete storagePtr_;
        };

    // Member Functions

        //- Storage
        virtual storage& data() const
        {
            if (!storagePtr_)
            {
                FatalErrorIn(type() + "::data()")
                    << "Storage pointer is not allocated. "
                        << "Most likely, this is not the master manager."
                        << abort(FatalError);
            }

            return *storagePtr_;
        };

        //- Solution control
        pimpleControl& pimple() const
        {
            return pimple_;
        };
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "interTrackManager.C"
#   include "interTrackManagerStorage.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
