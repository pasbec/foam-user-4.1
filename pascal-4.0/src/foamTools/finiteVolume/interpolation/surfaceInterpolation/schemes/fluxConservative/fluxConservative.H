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

Class
    Foam::fluxConservative

Description
    Flux-conserving interpolation scheme class, which interpolates sigma
    such that the face flux based on its gradient grad(phi), the corresponding
    diffusivity gamma and optionally an additionally given jump flux F will
    be preserved according to:

    (gamma_f*grad(phi)_f) * n_f + F/|S_f| = const.


SourceFiles
    fluxConservative.C

\*---------------------------------------------------------------------------*/

#ifndef fluxConservative_H
#define fluxConservative_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "linear.H"

// TODO: The diffusivity gamma is regarded as discontinous at every face!

// TODO: Tensorial diffusivity?

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fluxConservative Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fluxConservative
:
    public surfaceInterpolationScheme<Type>
{

    // Private Member Functions

        //- Name of the (scalar) diffusivity
        const volScalarField& gamma_;

        //- Jump flux field
        const GeometricField<Type, fvsPatchField, surfaceMesh>* jumpFluxPtr_;

        //- Disallow default bitwise assignment
        void operator=(const fluxConservative&);


public:

    //- Runtime type information
    TypeName("fluxConservative");


    // Constructors

        //- Construct from mesh and Istream
        fluxConservative
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            gamma_
            (
                mesh.thisDb().lookupObject<volScalarField>
                (
                    word(is)
                )
            ),
            jumpFluxPtr_(NULL)
        {
            if (!is.eof())
            {
                const GeometricField<Type, fvsPatchField, surfaceMesh>& jumpFlux =
                    mesh.thisDb().lookupObject
                    <
                        GeometricField<Type, fvsPatchField, surfaceMesh>
                    >
                    (
                        word(is)
                    );

                jumpFluxPtr_ = &jumpFlux;
            }
        }


        //- Construct from mesh, faceFlux and Istream
        fluxConservative
        (
            const fvMesh& mesh,
            const surfaceScalarField&,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            gamma_
            (
                mesh.thisDb().lookupObject<volScalarField>
                (
                    word(is)
                )
            ),
            jumpFluxPtr_(NULL)
        {
            if (!is.eof())
            {
                const GeometricField<Type, fvsPatchField, surfaceMesh>& jumpFlux =
                    mesh.thisDb().lookupObject
                    <
                        GeometricField<Type, fvsPatchField, surfaceMesh>
                    >
                    (
                        word(is)
                    );

                jumpFluxPtr_ = &jumpFlux;
            }
        }


        //- Construct from mesh, gamma and jump flux pointer
        fluxConservative
        (
            const fvMesh& mesh,
            const volScalarField& gamma,
            const GeometricField<Type, fvsPatchField, surfaceMesh>* jumpFluxPtr = NULL
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            gamma_(gamma),
            jumpFluxPtr_(jumpFluxPtr)
        {}


    // Member Functions

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return jumpFluxPtr_;
        }

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;

        //- Return the explicit correction to the face-interpolate
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fluxConservativeTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
