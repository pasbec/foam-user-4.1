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
    solverControl

Description
    ...

SourceFiles
    solverControl.C

\*---------------------------------------------------------------------------*/

#ifndef solverControl_H
#define solverControl_H

#include "foamTime.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define setControlScope(name)                                                \
                                                                              \
    name##Control& name = control.name();                                     \
    name##Control& control = name;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class solverControl Declaration
\*---------------------------------------------------------------------------*/

template <class MESH>
class solverControl
:
    public IOdictionary
{
protected:

    // Protected Data

        // Basic

            //- Reference to time
            const Time& time_;

            //- Reference to region mesh
            const MESH& mesh_;

            //- Name of base region
            const word baseRegionName_;

            //- Label of base region
            const label baseRegion_;

            //- Properties dictionary
            IOdictionary propDict_;

        // Dictionaries

            //- Solution dictionary (of base region)
            IOdictionary& baseSolutionDict_;

    // Protected Member Functions

private:

    // Private Data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        solverControl(const solverControl&);

        //- Disallow default bitwise assignment
        void operator=(const solverControl&);

public:

    // Constructors

        //- Construct from mesh and name
        solverControl
        (
            const MESH& mesh,
            const word& name = "solver",
            const label& base = -1
        );

    // Destructor

        virtual ~solverControl() {};

    // Member Functions

        // Basic

            //- Reference to mesh
            const Time& time() const
            {
                return time_;
            };

            //- Reference to mesh
            const MESH& mesh() const
            {
                return mesh_;
            };

            //- Reference properties dictionary
            const IOdictionary& properties() const
            {
                return propDict_;
            };

        // Properties

            //- Reference to base region name
            const word& baseName() const
            {
                return baseRegionName_;
            };

            //- Reference to base region ID
            const label& base() const
            {
                return baseRegion_;
            };

        // Time and database

            //- WriteData member function required by regIOobject
            bool writeData(Ostream&) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "solverControl.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
