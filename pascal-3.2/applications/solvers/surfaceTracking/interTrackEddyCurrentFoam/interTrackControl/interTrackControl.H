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

\*---------------------------------------------------------------------------*/

#ifndef interTrackControl_H
#define interTrackControl_H

#include "PstreamReduceOps.H"
#include "mathematicalConstants.H"
#include "physicalConstants.H"
#include "solverControl.H"
#include "regionModelling.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class interTrackControl Declaration
\*---------------------------------------------------------------------------*/

class interTrackControl
:
    public solverControl<regionFvMesh>
{
private:

    // Private Data

        // Controls

        // Properties

            //- Name of fluidA region
            word fluidAregionName_;

            //- Label of fluidA region
            label fluidAregion_;

        // Solution properties

        // Solution state variables

        // Solution residual variables

    // Private Member Functions

        //- Disallow default bitwise copy construct
        interTrackControl(const interTrackControl&);

        //- Disallow default bitwise assignment
        void operator=(const interTrackControl&);

        //- Inline helper

        //- Helper

public:

    //- Runtime type information
    TypeName("interTrackControl");

    // Constructors

        //- Construct from mesh and name
        interTrackControl
        (
            const regionFvMesh& rmesh,
            const word& name = "interTrack"
        );

    // Destructor

        virtual ~interTrackControl() {};

    // Member Functions

        // Properties

            //- Reference to fluidA region name
            const word& fluidAname() const
            {
                return fluidAregionName_;
            };

            //- Reference to fluidA region ID
            const label& fluidA() const
            {
                return fluidAregion_;
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
#   include "interTrackControl.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
