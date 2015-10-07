
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
    Foam::fluxConservativeValue

Description
    TODO

SourceFiles
    fluxConservativeValue.C

\*---------------------------------------------------------------------------*/

#ifndef fluxConservativeValue_H
#define fluxConservativeValue_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fluxConservativeValue Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fluxConservativeValue
:
    public surfaceInterpolationScheme<Type>
{
    // Private member data

        //- Conductivity reference
        const volScalarField& conductivity_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const fluxConservativeValue&);


public:

    //- Runtime type information
    TypeName("fluxConservativeValue");


    // Constructors

        //- Construct from conductivity
        fluxConservativeValue(const volScalarField& conductivity)
        :
            surfaceInterpolationScheme<Type>(conductivity.mesh()),
            conductivity_(conductivity)
        {
            checkConductivity();
        }

        //- Construct from Istream
        fluxConservativeValue(const fvMesh& mesh, Istream& is)
        :
            surfaceInterpolationScheme<Type>(mesh),
            conductivity_
            (
                mesh.lookupObject<volScalarField>
                (
                    word(is)
                )
            )
        {
            checkConductivity();
        }

        //- Construct from faceFlux and Istream
        fluxConservativeValue
        (
            const fvMesh& mesh,
            const surfaceScalarField&,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            conductivity_
            (
                mesh.lookupObject<volScalarField>
                (
                    word(is)
                )
            )
        {
            checkConductivity();
        }


    // Member Functions

        //- Check conductivity
        void checkConductivity() const
        {
            if (min(conductivity_).value() < 0.0)
            {
                FatalErrorIn
                (
                    "tmp<surfaceScalarField> fluxConservativeValue::weights\n"
                    "(\n"
                    "    const GeometricField<Type, fvPatchField, volMesh>&\n"
                    ")\n"
                )   << "Conductivity " << conductivity_.name()
                    << " is locally negative!" << endl
                    << "This probably happend due to an unbounded solution or."
                    << abort(FatalError);
            }
        }


        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
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
                    mesh_.weights()
                )
            );

            surfaceScalarField& w = tWeightingFactors();
            scalarField& wIn = w.internalField();

            const volScalarField& c = conductivity_;
            const scalarField& cIn = c.internalField();

            const unallocLabelList& own = mesh_.owner();
            const unallocLabelList& nei = mesh_.neighbour();

            const scalarField& wl = mesh_.weights();
            const scalarField& d = mesh_.deltaCoeffs();

            forAll (own, facei)
            {
                // Face centre to cell centre distances
                scalar fN = wl[facei] / d[facei];
                scalar Pf = 1.0 / d[facei] - fN;

                // New owner/neighbour coeffs
                scalar cOwn = cIn[own[facei]] / Pf;
                scalar cNei = cIn[nei[facei]] / fN;

                // TODO: What's wrong with this?
                wIn[facei] = cOwn/(cOwn + cNei + VSMALL);
            }

            return tWeightingFactors;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
