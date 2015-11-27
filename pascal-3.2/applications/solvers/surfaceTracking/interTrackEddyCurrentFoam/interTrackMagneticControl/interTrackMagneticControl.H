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
    interTrackMagneticControl

Description
    ...

SourceFiles
    interTrackMagneticControl.C

\*---------------------------------------------------------------------------*/

#ifndef interTrackMagneticControl_H
#define interTrackMagneticControl_H

#include "PstreamReduceOps.H"
#include "mathematicalConstants.H"
#include "physicalConstants.H"
#include "solverControl.H"
#include "regionModelling.H"

#include "interTrackControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class interTrackMagneticControl Declaration
\*---------------------------------------------------------------------------*/

class interTrackMagneticControl
:
    public solverControl<regionDynamicFvMesh>
{
private:

    // Private Data

        // Controls

            //- Interface-tracking specific controls
            mutable interTrackControl track_;

        // Properties

            //- Name of dynamic region
            word dynamicRegionName_;

            //- Label of dynamic region
            label dynamicRegion_;

        // Solution properties

        // Solution state variables

        // Solution residual variables

    // Private Member Functions

        //- Disallow default bitwise copy construct
        interTrackMagneticControl(const interTrackMagneticControl&);

        //- Disallow default bitwise assignment
        void operator=(const interTrackMagneticControl&);

        //- Inline helper

        //- Helper

public:

    //- Runtime type information
    TypeName("interTrackMagneticControl");

    // Constructors

        //- Construct from mesh and name
        interTrackMagneticControl
        (
            const regionDynamicFvMesh& rmesh,
            const word& name = "interTrackMagnetic"
        );

    // Destructor

        virtual ~interTrackMagneticControl() {};

    // Member Functions

        // Controls

            //- Reference to interface-tracking control
            interTrackControl& interTrack() const
            {
                return track_;
            };

        // Properties

            //- Reference to fluidA region name
            const word& fluidAname() const
            {
                return track_.fluidAname();
            };

            //- Reference to fluidA region ID
            const label& fluidA() const
            {
                return track_.fluidA();
            };

            //- Reference to fluidA region name
            const word& dynamicName() const
            {
                return dynamicRegionName_;
            };

            //- Reference to fluidA region ID
            const label& dynamic() const
            {
                return dynamicRegion_;
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
#   include "interTrackMagneticControl.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
