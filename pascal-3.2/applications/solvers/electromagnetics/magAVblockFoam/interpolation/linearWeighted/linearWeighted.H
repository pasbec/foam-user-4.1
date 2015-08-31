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

Class
    Foam::linearWeighted

Description
    TODO

SourceFiles
    linearWeighted.C

\*---------------------------------------------------------------------------*/

#ifndef linearWeighted_H
#define linearWeighted_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class linearWeighted Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class linearWeighted
:
    public surfaceInterpolationScheme<Type>
{
    // Private member data

        const volScalarField& weightWeightField_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const linearWeighted&);


public:

    //- Runtime type information
    TypeName("linearWeighted");


    // Constructors

        //- Construct from weightWeightField
        linearWeighted(const volScalarField& weights)
        :
            surfaceInterpolationScheme<Type>(weights.mesh()),
            weightWeightField_(weights)
        {}

        //- Construct from Istream
        linearWeighted(const fvMesh& mesh, Istream& is)
        :
            surfaceInterpolationScheme<Type>(mesh),
            weightWeightField_
            (
                this->mesh().objectRegistry::template
                lookupObject<const volScalarField>(word(is))
            )
        {}

        //- Construct from faceFlux and Istream
        linearWeighted
        (
            const fvMesh& mesh,
            const surfaceScalarField&,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            weightWeightField_
            (
                this->mesh().objectRegistry::template
                lookupObject<const volScalarField>
                (
                    word(is)
                )
            )
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            // TODO: mesh_ of parent class does not work?
            const fvMesh& mesh_ = this->mesh();

            // Copy linear weighting factors
            tmp<surfaceScalarField> tWeightingFactors
            (
                new surfaceScalarField
                (
                    IOobject
                    (
                        "weightingFactors",
                        mesh_.pointsInstance(),
                        mesh_
                    ),
                    mesh_.surfaceInterpolation::weights()
                )
            );

            // Short refs
            surfaceScalarField& weightingFactors = tWeightingFactors();
            scalarField& w = weightingFactors.internalField();

            // Short ref to weightWeightField
            const volScalarField& wwf = weightWeightField_;
            const scalarField& wwfi = wwf.internalField();

            // Set local references to mesh data
            const unallocLabelList& owner = mesh_.owner();
            const unallocLabelList& neighbour = mesh_.neighbour();

            const vectorField& Cf = mesh_.faceCentres();
            const vectorField& C = mesh_.cellCentres();
            const vectorField& Sf = mesh_.faceAreas();

            forAll (owner, facei)
            {
                scalar SfdOwn = mag(Sf[facei] & (Cf[facei] - C[owner[facei]]));
                scalar SfdNei = mag(Sf[facei] & (C[neighbour[facei]] - Cf[facei]));

                scalar fOwn = wwfi[owner[facei]];
                scalar fNei = wwfi[neighbour[facei]];

                w[facei] = fOwn*SfdNei/(fNei*SfdOwn + fOwn*SfdNei);
//                 w[facei] = SfdNei/fNei / (SfdOwn/fOwn + SfdNei/fNei);
            }

            return tWeightingFactors;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
