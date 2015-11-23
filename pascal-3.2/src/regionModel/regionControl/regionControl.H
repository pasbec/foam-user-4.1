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
    regionControl

Description
    ...

SourceFiles
    regionControl.C

\*---------------------------------------------------------------------------*/

#ifndef regionControl_H
#define regionControl_H

#include "foamTime.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class regionControl Declaration
\*---------------------------------------------------------------------------*/

template <class RegionGeoMesh>
class regionControl
:
    regIOobject
{
public:

    // Public typedefs

        typedef typename RegionGeoMesh::Mesh RegionMesh;

protected:

    // Protected Data

    // Protected Member Functions

private:

    // Private Data

        //- Reference to region mesh
        const RegionMesh& rMesh_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionControl(const regionControl&);

        //- Disallow default bitwise assignment
        void operator=(const regionControl&);

public:

    //- Runtime type information
    TypeName("regionControl");

    // Constructors

        //- Construct from region mesh
        regionControl(const RegionMesh& rmesh);

    // Destructor

        virtual ~regionControl() {};

    // Member Functions

        // Access

            //- Reference to mesh
            const RegionMesh& rmesh() const
            {
                return rMesh_;
            };

            //- WriteData member function required by regIOobject
            bool writeData(Ostream&) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "regionControl.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
