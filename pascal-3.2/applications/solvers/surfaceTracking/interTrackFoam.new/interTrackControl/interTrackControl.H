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
    interTrackControl

Description
    ...

SourceFiles
    interTrackControl.C
    interTrackControlStorage.C

\*---------------------------------------------------------------------------*/

#ifndef interTrackControl_H
#define interTrackControl_H

#include "solverControl.H"

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "trackedSurface.H"

// #include "PstreamReduceOps.H"
// #include "mathematicalConstants.H"
// #include "physicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class interTrackControl Declaration
\*---------------------------------------------------------------------------*/

class interTrackControl
:
    public solverControl<dynamicFvMesh>
{
public:

    // Public Nested Classes

        //- Storage
        struct storage
        :
            solverControl<dynamicFvMesh>::storage
        {
        private:

            // Private Data

                //- Reference to control
                interTrackControl& control_;

            // Private Member Functions

                //- Disallow default bitwise copy construct
                storage(const storage&);

                //- Disallow default bitwise assignment
                void operator=(const storage&);

                //- Return reference to control
                interTrackControl& control()
                {
                    return control_;
                }

                //- Initialize data sequence
                void init();

            // Macro-generated data entries with accessors

                //- Pressure
                solverControlStorage_add(volScalarField, p);

                //- Velocity
                solverControlStorage_add(volVectorField, U);

                //- Face flux
                solverControlStorage_add(surfaceScalarField, phi);

                //- Density
                solverControlStorage_add(volScalarField, rho);

                //- Dynamic viscosity
                solverControlStorage_add(volScalarField, mu);

                //- Fluid indicator
                solverControlStorage_add(volScalarField, fluidIndicator);

                //- Tracked surface
                solverControlStorage_add(trackedSurface, interface);

        public:

            // Constructor

                //- Construct from control
                storage(interTrackControl& control);
        };

private:

    // Private Data

        //- Storage
        storage* storagePtr_;

        //- Controls

        //- Properties

        //- Dictionaries

        //- Solution properties

        //- Solution state variables

        //- Solution residual variables

    // Private Member Functions

        //- Disallow default bitwise copy construct
        interTrackControl(const interTrackControl&);

        //- Disallow default bitwise assignment
        void operator=(const interTrackControl&);

        //- Inline helper

        //- Helper

public:

    // Runtime type information
    TypeName("interTrackControl");

    // Constructors

        //- Construct from mesh and name
        interTrackControl
        (
            const argList& args,
            Time& time,
            dynamicFvMesh& mesh,
            const word& name = "interTrack",
            const bool& master = true
        );

    // Destructor

        virtual ~interTrackControl()
        {
            if (master()) delete storagePtr_;
        };

    // Member Functions

        //- Storage
        virtual storage& data()
        {
            return *storagePtr_;
        };

        //- Controls

        //- Properties

        //- Dictionaries

        //- Solution properties

        //- Solution state variables

        //- Solution residual variables

        //- Time and database

            //- WriteData member function required by regIOobject
            bool writeData(Ostream&) const
            {
                return false;
            };
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "interTrackControl.C"
#   include "interTrackControlStorage.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
