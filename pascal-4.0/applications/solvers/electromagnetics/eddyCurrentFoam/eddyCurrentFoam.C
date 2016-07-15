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
    eddyCurrentFoam

Description
    ...

\*---------------------------------------------------------------------------*/

#include "eddyCurrentApp.H"

// TODO: Derived gradient boundary condition for VRe/VIm in conductor region
// TODO: Remove interfaceLabel? Use alias of calculatedFvPatchField or special patchField

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    using namespace Foam;

    argList::validOptions.insert("overwrite", "");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createRegionMeshUninitialized.H"

    using namespace eddyCurrentApp;
    using namespace eddyCurrentApp::Region;

    Manager manager(args, runTime, regionMesh);

    SM_GLOBALREGIONSCOPE(DEFAULT);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    uniformDimensionedScalarField& omega0 = globalStorage.omega0();

    while(manager.once())
    {
#       include "A0BiotSavart.H"

#       include "AVLoop.H"

#       include "BUpdate.H"

#       include "derivedFields.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    return(0);
}

// ************************************************************************* //
