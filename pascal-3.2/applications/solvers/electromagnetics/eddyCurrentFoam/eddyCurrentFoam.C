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

Application
    ...

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvBlockMatrix.H"
#include "eddyCurrentControl.H"

// TODO: Derived gradient boundary condition for VRe/VIm in conductor region

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validOptions.insert("overwrite", "");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createRegionMesh.H"

#   include "eddyCurrentRegionFields.H"
#   include "eddyCurrentBaseFields.H"
#   include "eddyCurrentConductorFields.H"

#ifndef namespaceFoam
#define namespaceFoam
    using namespace Foam;
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    if (!args.optionFound("overwrite"))
    {
        runTime++;

        control.msg().time();
    }

#   include "eddyCurrentAVloop.H"

#   include "eddyCurrentDerivedFields.H"

    runTime.writeNow();

    control.msg().executionTime();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    control.msg().end();

    return(0);
}

// ************************************************************************* //
