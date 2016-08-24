/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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
    makeRegionToRegionAddressing

Description
    Make (calculate if not already present) region-to-region addressings
    for multi-region purposes.

Usage
    - makeRegionToRegionAddressing -fromRegions wordList -toRegions wordList \n
    Pre-calculate region-to-region addressings for multi-region purposes

    @param -fromRegions wordList \n
    Specify particular regions as source regions

    @param -toRegions wordList \n
    Specify particular regions as target regions

    E.g.
    - makeRegionToRegionAddressing -fromRegions "(region0 region_fluid)"
                                   -toRegions "(region0 region_fluid)"


Author
    Pascal Beckstein

// TODO: Force presence of parameters

// TODO: Check sense and dependency of from and to

\*---------------------------------------------------------------------------*/

#include "parRun.H"

#include "objectRegistry.H"
#include "foamTime.H"
#include "argList.H"

#include "OSspecific.H"
#include "argList.H"
#include "timeSelector.H"

#include "addressingTypes.H"
#include "regionModelling.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

int main(int argc, char *argv[])
{
    // Add list for source regions
    argList::validOptions.insert("fromRegions", "wordList");

    // Add list for target regions
    argList::validOptions.insert("toRegions", "wordList");

    #include "setRootCase.H"
    #include "createTime.H"

    #include "createRegionPolyMeshUninitialized.H"

    wordList fromRegionNames(args.optionLookup("fromRegions")());
    wordList toRegionNames(args.optionLookup("toRegions")());

    regionMesh.init(fromRegionNames);

    forAll (regionMesh, fromRegionI)
    {
        forAll (toRegionNames, toRegionI)
        {
            // Loop over point, face, cell, boundary
            wordList addressingTypeNames =
                addressingTypes::addressingTypeName.toc();

            // Calc addressings
            forAll (addressingTypeNames, addressingTypeI)
            {
                word addressingTypeName = addressingTypeNames[addressingTypeI];

                regionMesh.typeAddressing
                (
                    addressingTypes::addressingTypeName[addressingTypeName],
                    regionMesh.regions()[fromRegionI],
                    toRegionNames[toRegionI]
                );
            }
        }
    }

    Info << "\nEnd\n" << endl;

    return(0);
}

// ************************************************************************* //

