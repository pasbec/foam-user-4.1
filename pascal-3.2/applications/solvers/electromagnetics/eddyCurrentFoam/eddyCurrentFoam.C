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

#include "eddyCurrentManager.H"

#ifndef namespaceFoam
#define namespaceFoam
    using namespace Foam;
#endif

// TODO: Derived gradient boundary condition for VRe/VIm in conductor region

// TODO: Extrapolation sometimes completely fails

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validOptions.insert("overwrite", "");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createRegionMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    eddyCurrentManager manager(args, runTime, regionMesh);

    SM_REGIONSCOPE_DEFAULTREGION(eddyCurrentManager);

    uniformDimensionedScalarField& omega0 = manager.storage().omega0();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    while(manager.once())
    {
// #       include "AVloop.H"

        {
//             typedef eddyCurrentManager Manager;
//             typedef Manager::Control GlobalControl;
//             typedef Manager::Settings GlobalSettings;
//             typedef Manager::Storage GlobalStorage;
// 
//             GlobalControl& globalControl = manager.control();
//             GlobalSettings& globalSettings = manager.settings();
//             GlobalStorage& globalStorage = manager.storage();
// 
//             //typedef typename Manager::Region<manager.regions().defaultRegion()>::Type Region;
//             typedef Manager::DefaultRegion Region;
//             typedef Region::Mesh Mesh;
//             typedef Region::Control Control;
//             typedef Region::Settings Settings;
//             typedef Region::Storage Storage;
// 
//             Region& region = manager.regions().defaultRegion();
//             Mesh& mesh = region.mesh();
//             Control& control = region.control();
//             Settings& settings = region.settings();
//             Storage& storage = region.storage();

//             volScalarField& sigma = manager.storage().sigma()[mesh.name()];
//
//             volVectorField& ARe = manager.storage().ARe()[mesh.name()];
//             volVectorField& AIm = manager.storage().AIm()[mesh.name()];
//
//             volVectorField& VReGrad = manager.storage().VReGrad()[mesh.name()];
//             volVectorField& VImGrad = manager.storage().VImGrad()[mesh.name()];
//
//             volVectorField& FL = manager.storage().FL()[mesh.name()];
//             volScalarField& pB = manager.storage().pB()[mesh.name()];
//
//             volVectorField& jRe = manager.regions().defaultRegion().storage().jRe();
//             volVectorField& jIm = manager.regions().defaultRegion().storage().jIm();
//             volVectorField& BRe = manager.regions().defaultRegion().storage().BRe();
//             volVectorField& BIm = manager.regions().defaultRegion().storage().BIm();
//
//             regionVolVectorField& FL = manager.storage().FL();
//             regionVolScalarField& pB = manager.storage().pB();

// #           include "derivedFields.H"
        }
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    return(0);
}

// ************************************************************************* //
