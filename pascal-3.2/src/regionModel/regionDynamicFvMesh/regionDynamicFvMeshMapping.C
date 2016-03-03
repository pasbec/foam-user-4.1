/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

\*---------------------------------------------------------------------------*/

#include "tetFemMatrices.H"
#include "tetPointFields.H"
#include "fixedValueTetPolyPatchFields.H"
#include "faceTetPolyPatch.H"
#include "tetPolyPatchInterpolation.H"
#include "fixedValuePointPatchFields.H"
#include "twoDPointCorrector.H"

#include "regionDynamicFvMesh.H"

// TODO [High]: What if directMapped patches have been manually modified be changeDictionary?

// TODO [High]: What happens on parallel ( hint: mappedpatchBase::distribute() )

// TODO [High]: What happens on parallel ( remember "preserveFaceZones" )

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void regionDynamicFvMesh::patchMapMeshVelocityDirectMapped
(
    const label& fromRegionI,
    const label& toRegionI
) const
{
    if (isFeMotionSolver(toRegionI))
    {
        tetPointVectorField& toMotionU =
        const_cast<tetPointVectorField&>
        (
            mesh(toRegionI).objectRegistry::
            lookupObject<tetPointVectorField>
            (
                "motionU"
            )
        );

        labelList patchMap =
            regionPolyMesh::patchMapDirectMapped
            (
                fromRegionI,
                toRegionI
            )[0];

        forAll(mesh(fromRegionI).boundary(), fromPatchI)
        {
            label toPatchI = patchMap[fromPatchI];

            if (toPatchI != -1)
            {
                const polyPatch& fromPolyPatch =
                    mesh(fromRegionI).boundaryMesh()[fromPatchI];

                const polyPatch& toPolyPatch =
                    mesh(toRegionI).boundaryMesh()[toPatchI];

                tmp<vectorField> ttoTotalDisplacement
                (
                    new vectorField
                    (
                        fromPolyPatch.localPoints()
                      - toPolyPatch.localPoints()
                    )
                );

                vectorField& toTotalDisplacement = ttoTotalDisplacement();

                fixedValueTetPolyPatchVectorField& toMotionUpatch =
                    refCast<fixedValueTetPolyPatchVectorField>
                    (
                        toMotionU.boundaryField()[toPatchI]
                    );

                tetPolyPatchInterpolation toTppiPatch
                (
                    refCast<const faceTetPolyPatch>
                    (
                        toMotionUpatch.patch()
                    )
                );

                toMotionUpatch ==
                    toTppiPatch.pointToPointInterpolate
                    (
                        toTotalDisplacement/time().deltaT().value()
                    );

                ttoTotalDisplacement.clear();
            }
        }
    }
    else if (isFvMotionSolver(toRegionI))
    {
        // Mesh velocity field
        pointVectorField& toPointMotionU =
        const_cast<pointVectorField&>
        (
            mesh(toRegionI).objectRegistry::
            lookupObject<pointVectorField>
            (
                "pointMotionU"
            )
        );

        labelList patchMap =
            patchMapDirectMapped
            (
                fromRegionI,
                toRegionI
            )[0];

        forAll(mesh(fromRegionI).boundary(), fromPatchI)
        {
            label toPatchI = patchMap[fromPatchI];

            if (toPatchI != -1)
            {
                const polyPatch& fromPolyPatch =
                    mesh(fromRegionI).boundaryMesh()[fromPatchI];

                const polyPatch& toPolyPatch =
                    mesh(toRegionI).boundaryMesh()[toPatchI];

                tmp<vectorField> ttoTotalDisplacement
                (
                    new vectorField
                    (
                        fromPolyPatch.localPoints()
                      - toPolyPatch.localPoints()
                    )
                );

                vectorField& toTotalDisplacement = ttoTotalDisplacement();

                fixedValuePointPatchVectorField& toPointMotionUpatch =
                    refCast<fixedValuePointPatchVectorField>
                    (
                        toPointMotionU.boundaryField()[toPatchI]
                    );

                toPointMotionUpatch ==
                    toTotalDisplacement/time().deltaT().value();

                ttoTotalDisplacement.clear();
            }
        }
    }
}

void regionDynamicFvMesh::patchMapMeshVelocityDirectMapped
(
    const word& fromRegionName,
    const word& toRegionName
) const
{
    label fromRegionI = regionIndex(fromRegionName);
    label toRegionI = regionIndex(toRegionName);

    return patchMapMeshVelocityDirectMapped
    (
        fromRegionI,
        toRegionI
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

