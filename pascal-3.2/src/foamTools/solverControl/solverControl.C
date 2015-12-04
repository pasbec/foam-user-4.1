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

\*---------------------------------------------------------------------------*/

#include "solverControl.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class MESH>
solverControl<MESH>::solverControl
(
    const MESH& mesh,
    const word& name,
    const bool& master,
    const label& regionI0
)
:
    IOdictionary
    (
        IOobject
        (
            name + "Control",
            mesh.time().timeName(),
            "uniform",
            mesh.time().db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    time_(mesh.time()),
    mesh_(mesh),
    master_(master),
    baseRegionName_(polyMesh::defaultRegion),
    baseRegion_(regionI0),
    propDict_
    (
        IOdictionary
        (
            IOobject
            (
                name + "Properties",
                time_.constant(),
                time_.db(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        )
    ),
    baseSolutionDict_(mesh_.solutionDict()),
    msg_(NULL)
{
    if (master) msg_ = new solverControlMessages(mesh_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class MESH>
bool solverControl<MESH>::writeData(Ostream&) const
{
    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

