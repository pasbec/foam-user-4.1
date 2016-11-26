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
    Foam::fv::jumpGaussLaplacian

Description
    ...

SourceFiles
    jumpGaussLaplacian.C

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#ifndef jumpGaussLaplacian_H
#define jumpGaussLaplacian_H

#include "laplacianScheme.H"
#include "harmonic.H"
#include "snGradScheme.H"
#include "uncorrectedSnGrad.H"

// TODO: The diffusivity gamma is regarded as discontinous at every face!

// TODO: Tensorial diffusivity?

// TODO: Think about CORRECT implementation with snGrad correction. The
//       corrected (explicit) gradient may not be interpolated from the cell
//       values using linear interpolation!!! We MUST use fluxConservative
//       or a corresponding interpolation scheme to calculate the gradient
//       correction in normal direction.

// TODO: Proper template specialisation? The diffusivity may be
//       any GType. Specialization is necessary for e.g. tensors. Currently
//       we simply use the mag-function as workaround and in order to avoid
//       having to create these.

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class jumpGaussLaplacian Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GType>
class jumpGaussLaplacian
:
    public fv::laplacianScheme<Type, GType>
{

private:

    // Private data

        //- Name of the (scalar) diffusivity
        const GeometricField<GType, fvPatchField, volMesh>& gamma_;

        //- Jump flux field
        const GeometricField<Type, fvsPatchField, surfaceMesh>* jumpFluxPtr_;

    // Private Member Functions

        tmp<fvMatrix<Type> > fvmLaplacianUncorrected
        (
            const surfaceScalarField& gammaMagSf,
            const word& vfGammaName,
            const GeometricField<Type, fvPatchField, volMesh>& vf
        );

        void fvmLaplacianJumpFlux
        (
            fvMatrix<Type>& fvm,
            const surfaceScalarField& gammaMagSf,
            const GeometricField<Type, fvPatchField, volMesh>& vf
        );

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > gammaSnGradCorr
        (
            const surfaceVectorField& SfGammaCorr,
            const word& vfGammaName,
            const GeometricField<Type, fvPatchField, volMesh>& vf
        );

        //- Disallow default bitwise copy construct
        jumpGaussLaplacian(const jumpGaussLaplacian&);

        //- Disallow default bitwise assignment
        void operator=(const jumpGaussLaplacian&);


public:

    //- Runtime type information
    TypeName("jumpGauss");

    // Constructors

        //- Construct from Istream
        jumpGaussLaplacian
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            laplacianScheme<Type, GType>
            (
                mesh,
                tmp<surfaceInterpolationScheme<GType> >
                (
                    new harmonic<GType>(mesh)
                ),
// TODO FIXME: Special snGradSchemes are necessary?
                tmp<snGradScheme<Type> >
                (
                    new uncorrectedSnGrad<Type>(mesh)
                )
            ),
            gamma_
            (
                mesh.thisDb().lookupObject
                <
                    GeometricField<GType, fvPatchField, volMesh>
                >
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

            word schemeSigma = mesh.schemesDict().interpolationScheme
            (
                "interpolate(" + gamma_.name() + ")"
            );

            if (schemeSigma != "harmonic")
            {
                FatalErrorIn
                (
                    "jumpGaussLaplacian::jumpGaussLaplacian(...)"
                )   << "For this Laplacian, gamma needs to be"
                    << " interpolated harmonically!"
                    << exit(FatalError);
            }
        }

    // Destructor

        virtual ~jumpGaussLaplacian()
        {}

    // Member Functions

        // -Laplacian method with name of volField of gamma
        tmp<fvMatrix<Type> > fvmLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            const word&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        // -Laplacian method with surfaceField of gamma
        virtual tmp<fvMatrix<Type> > fvmLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        // -Laplacian method which interpolates gamma and passes its name
        virtual tmp<fvMatrix<Type> > fvmLaplacian
        (
            const GeometricField<GType, fvPatchField, volMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        // -NOTE: Deactivated!
        virtual tmp<GeometricField<Type, fvPatchField, volMesh> > fvcLaplacian
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        // -Laplacian calculation with name of volField of gamma
        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            const word&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        // -Laplacian calculation with surfaceField of gamma
        virtual tmp<GeometricField<Type, fvPatchField, volMesh> > fvcLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        // -Laplacian calculation which interpolates gamma and passes its name
        virtual tmp<GeometricField<Type, fvPatchField, volMesh> > fvcLaplacian
        (
            const GeometricField<GType, fvPatchField, volMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );
};


// Use macros to emulate partial-specialisation of the the Laplacian functions
// for scalar diffusivity gamma

#define defineLaplacianScalarGamma(Type)                                    \
                                                                            \
template<>                                                                  \
tmp<fvMatrix<Type> >                                                        \
jumpGaussLaplacian<Type, scalar>::fvmLaplacian                              \
(                                                                           \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>&,              \
    const word&,                                                            \
    const GeometricField<Type, fvPatchField, volMesh>&                      \
);                                                                          \
                                                                            \
template<>                                                                  \
tmp<fvMatrix<Type> >                                                        \
jumpGaussLaplacian<Type, scalar>::fvmLaplacian                              \
(                                                                           \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>&,              \
    const GeometricField<Type, fvPatchField, volMesh>&                      \
);                                                                          \
                                                                            \
template<>                                                                  \
tmp<GeometricField<Type, fvPatchField, volMesh> >                           \
jumpGaussLaplacian<Type, scalar>::fvcLaplacian                              \
(                                                                           \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>&,              \
    const word&,                                                            \
    const GeometricField<Type, fvPatchField, volMesh>&                      \
);                                                                          \
                                                                            \
template<>                                                                  \
tmp<GeometricField<Type, fvPatchField, volMesh> >                           \
jumpGaussLaplacian<Type, scalar>::fvcLaplacian                              \
(                                                                           \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>&,              \
    const GeometricField<Type, fvPatchField, volMesh>&                      \
);


defineLaplacianScalarGamma(scalar);
defineLaplacianScalarGamma(vector);
defineLaplacianScalarGamma(sphericalTensor);
defineLaplacianScalarGamma(symmTensor);
defineLaplacianScalarGamma(tensor);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "jumpGaussLaplacian.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
