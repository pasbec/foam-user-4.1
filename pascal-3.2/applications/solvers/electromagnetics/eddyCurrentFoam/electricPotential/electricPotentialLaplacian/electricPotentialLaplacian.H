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
    Foam::fv::electricPotentialLaplacian

Description
    ...

SourceFiles
    electricPotentialLaplacian.C

// TODO: Currently it is NOT SUPPORTED or even tested to put the
//       conducting region on any type of boundary

\*---------------------------------------------------------------------------*/

#ifndef electricPotentialLaplacian_H
#define electricPotentialLaplacian_H

#include "laplacianScheme.H"
#include "electricPotentialLaplacianBase.H"
#include "harmonic.H"
#include "snGradScheme.H"
#include "uncorrectedSnGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class electricPotentialLaplacian Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GType>
class electricPotentialLaplacian
:
    public fv::electricPotentialLaplacianBase,
    public fv::laplacianScheme<Type, GType>
{

private:

    // Private data

        //- Name of the field of rate of change of
        //  flux of magnetic vector potential
        word nameDdtAflux_;

        //- Name of the faceSet for jumps of conductivity
        word nameInterface_;

    // Private Member Functions

        tmp<fvMatrix<Type> > fvmLaplacianUncorrected
        (
            const surfaceScalarField& gammaMagSf,
            const word&,
            GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<Field<Type> > fvmLaplacianSource
        (
            const surfaceScalarField& gammaMagSf,
            const word&
        );

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > gammaSnGradCorr
        (
            const surfaceVectorField& SfGammaCorr,
            const word&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        //- Disallow default bitwise copy construct
        electricPotentialLaplacian(const electricPotentialLaplacian&);

        //- Disallow default bitwise assignment
        void operator=(const electricPotentialLaplacian&);


public:

    //- Runtime type information
    TypeName("electricPotentialLaplacian");

    // Constructors

        //- Construct null
        electricPotentialLaplacian
        (
            const fvMesh& mesh
        )
        :
            electricPotentialLaplacianBase(GAUSS),
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
            nameDdtAflux_("ddtAflux"),
            nameInterface_()
        {
            checkBaseScheme();
        }

        //- Construct from Istream
        electricPotentialLaplacian
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            electricPotentialLaplacianBase(is),
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
            nameDdtAflux_(is),
            nameInterface_()
        {
            checkBaseScheme();

            // TODO: Apply only for given faceSet?
            if (!is.eof())
            {
                nameInterface_ = word(is);
            }
        }

    // Destructor

        virtual ~electricPotentialLaplacian()
        {}

    // Member Functions

        // -Laplacian method with name of volField of gamma
        tmp<fvMatrix<Type> > fvmLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            const word&,
            GeometricField<Type, fvPatchField, volMesh>&
        );

        // -NOTE: Deactivated!
        virtual tmp<fvMatrix<Type> > fvmLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            GeometricField<Type, fvPatchField, volMesh>&
        );

        // -Laplacian method which interpolates gamma and passes its name
        virtual tmp<fvMatrix<Type> > fvmLaplacian
        (
            const GeometricField<GType, fvPatchField, volMesh>&,
            GeometricField<Type, fvPatchField, volMesh>&
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

        // -NOTE: Deactivated!
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
electricPotentialLaplacian<Type, scalar>::fvmLaplacian                      \
(                                                                           \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>&,              \
    const word&,                                                            \
    GeometricField<Type, fvPatchField, volMesh>&                            \
);                                                                          \
                                                                            \
template<>                                                                  \
tmp<fvMatrix<Type> >                                                        \
electricPotentialLaplacian<Type, scalar>::fvmLaplacian                      \
(                                                                           \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>&,              \
    GeometricField<Type, fvPatchField, volMesh>&                            \
);                                                                          \
                                                                            \
template<>                                                                  \
tmp<GeometricField<Type, fvPatchField, volMesh> >                           \
electricPotentialLaplacian<Type, scalar>::fvcLaplacian                      \
(                                                                           \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>&,              \
    const word&,                                                            \
    const GeometricField<Type, fvPatchField, volMesh>&                      \
);                                                                          \
                                                                            \
template<>                                                                  \
tmp<GeometricField<Type, fvPatchField, volMesh> >                           \
electricPotentialLaplacian<Type, scalar>::fvcLaplacian                      \
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
#   include "electricPotentialLaplacian.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
