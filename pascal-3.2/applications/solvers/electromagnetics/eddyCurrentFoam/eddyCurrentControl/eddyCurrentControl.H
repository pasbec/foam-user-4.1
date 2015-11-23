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
    eddyCurrentControl

Description
    ...

SourceFiles
    eddyCurrentControl.C

\*---------------------------------------------------------------------------*/

#ifndef eddyCurrentControl_H
#define eddyCurrentControl_H

#include "regionControls.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class eddyCurrentControl Declaration
\*---------------------------------------------------------------------------*/

class eddyCurrentControl
:
    regIOobject
{
private:

    // Private Data

        // Reference to region control
        const regionFvControl& control_;

        // Reference to region mesh
        const regionFvMesh& rmesh_;

        // Properties

            //- Properties dictionary
            IOdictionary propDict_;

            //- Name of base region
            word baseRegionName_;

            //- Name of conducting region
            word conductorRegionName_;

            //- Label of base region
            label baseRegion_;

            //- Label of conducting region
            label conductorRegion_;

            //- Frequency
            dimensionedScalar frequency_;

            //- Cyclic velocity based on frequency
            dimensionedScalar omega_;

        // Solution

            //- AV-solution dictionary
            dictionary& AVdict_;

            //- A-solution dictionary
            dictionary& Adict_;

            //- V-solution dictionary
            dictionary& Vdict_;

            //- Minimum number of iterations
            int minIter_;

            //- Maximum number of iterations
            int maxIter_;

            //- Solution tolerance
            scalar tol_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        eddyCurrentControl(const eddyCurrentControl&);

        //- Disallow default bitwise assignment
        void operator=(const eddyCurrentControl&);

public:

    //- Runtime type information
    TypeName("eddyCurrentControl");

    // Constructors

        //- Construct from regionFvControl
        eddyCurrentControl(const regionFvControl& rfvc);

    // Destructor

        virtual ~eddyCurrentControl() {};

    // Member Functions

        // Properties

            //- Reference to base ID
            const dictionary& properties() const
            {
                return propDict_;
            };

            //- Reference to base ID
            const label& baseRegion() const
            {
                return baseRegion_;
            };

            //- Reference to base ID
            const label& conductorRegion() const
            {
                return conductorRegion_;
            };

            //- Reference to region control
            const dimensionedScalar& frequency() const
            {
                return frequency_;
            };

            //- Reference to region control
            const dimensionedScalar& omega() const
            {
                return omega_;
            };

        // Solution

            //- Reference to AV-solution dictionary
            const dictionary& AVdict() const
            {
                return AVdict_;
            };

            //- Reference to A-solution dictionary
            const dictionary& Adict() const
            {
                return Adict_;
            };

            //- Reference to V-solution dictionary
            const dictionary& Vdict() const
            {
                return Vdict_;
            };

            //- Reference to minimum number of iterations
            const int& minIter() const
            {
                return minIter_;
            };

            //- Reference to maximum number of iterations
            const int& maxIter() const
            {
                return maxIter_;
            };

            //- Reference to solution tolerance
            const scalar& tol() const
            {
                return tol_;
            };

        // Other

            //- WriteData member function required by regIOobject
            bool writeData(Ostream&) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
