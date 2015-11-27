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
    interTrackEddyCurrentControl

Description
    ...

SourceFiles
    interTrackEddyCurrentControl.C

\*---------------------------------------------------------------------------*/

#ifndef interTrackEddyCurrentControl_H
#define interTrackEddyCurrentControl_H

#include "PstreamReduceOps.H"
#include "mathematicalConstants.H"
#include "physicalConstants.H"
#include "solverControl.H"
#include "regionModelling.H"

#include "interTrackMagneticControl.H"
#include "eddyCurrentControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class interTrackEddyCurrentControl Declaration
\*---------------------------------------------------------------------------*/

class interTrackEddyCurrentControl
:
    public solverControl<regionDynamicFvMesh>
{
private:

    // Private Data

        // Controls

            //- TODO: Time controls

            //- Interface-tracking specific controls
            mutable interTrackMagneticControl tramag_;

            //- Eddy current specific controls
            mutable eddyCurrentControl eddy_;

        // Properties

        // Solution properties

        // Solution state variables

        // Solution residual variables

    // Private Member Functions

        //- Disallow default bitwise copy construct
        interTrackEddyCurrentControl(const interTrackEddyCurrentControl&);

        //- Disallow default bitwise assignment
        void operator=(const interTrackEddyCurrentControl&);

        //- Inline helper

        //- Helper

public:

    //- Runtime type information
    TypeName("interTrackEddyCurrentControl");

    // Constructors

        //- Construct from mesh and name
        interTrackEddyCurrentControl
        (
            const regionDynamicFvMesh& rmesh,
            const word& name = "interTrackEddyCurrent"
        );

    // Destructor

        virtual ~interTrackEddyCurrentControl() {};

    // Member Functions

        // Controls

            //- Reference to interface-tracking control
            interTrackMagneticControl& interTrackMagnetic() const
            {
                return tramag_;
            };

            //- Reference to eddy current control
            eddyCurrentControl& eddyCurrent() const
            {
                return eddy_;
            };

        // Properties

            //- Reference to fluidA region name
            const word& fluidAname() const
            {
                return tramag_.fluidAname();
            };

            //- Reference to fluidA region ID
            const label& fluidA() const
            {
                return tramag_.fluidA();
            };

            //- Reference to dynamic region name
            const word& dynamicName() const
            {
                return tramag_.dynamicName();
            };

            //- Reference to dynamic region ID
            const label& dynamic() const
            {
                return tramag_.dynamic();
            };

            //- Reference to conductor region name
            const word& conductorName() const
            {
                return eddy_.conductorName();
            };

            //- Reference to conductor region ID
            const label& conductor() const
            {
                return eddy_.conductor();
            };

            //- Reference to region control
            const dimensionedScalar& frequency() const
            {
                return eddy_.frequency();
            };

            //- Reference to region control
            const dimensionedScalar& omega() const
            {
                return eddy_.omega();
            };

        // Solution properties

        // Solution state variables

        // Solution residual variables

        // Time and database

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
#   include "interTrackEddyCurrentControl.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
