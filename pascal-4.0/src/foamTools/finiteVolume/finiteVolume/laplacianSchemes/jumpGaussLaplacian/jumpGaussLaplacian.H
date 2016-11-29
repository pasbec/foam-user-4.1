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
#include "fluxConservative.H"
#include "harmonic.H"
#include "snGradScheme.H"
#include "uncorrectedSnGrad.H"

// TODO: Two-sided snGradSchemes necessary!

// TODO: Only include headers which are really needed!

// TODO: The diffusivity gamma is currently regarded as
//       discontinous at every face!

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

// TODO: Two-sided snGradSchemes necessary!
        //- Which snGradScheme is meant to be used. Just a dummy to do
        //  correction here
        word snGradSchemeName_;

// TODO: Two-sided snGradSchemes necessary!
        //- Is this scheme meant to be corrected?
        bool corrected_;

// TODO: Two-sided snGradSchemes necessary!
        //- Is this scheme's correction meant to be limited?
        scalar limiter_;

        //- Name of the (scalar) diffusivity
        const GeometricField<GType, fvPatchField, volMesh>& gamma_;

        //- Jump flux field
        const GeometricField<Type, fvsPatchField, surfaceMesh>* jumpFluxPtr_;

        //- Flux-conservative interpolation scheme
        const fluxConservative<Type> vfInterpolScheme_;

    // Private Member Functions

        //- Check for valid schmeme settings
        void checkSchemes() const;

        //- Uncorrected Laplacian
        tmp<fvMatrix<Type> > fvmLaplacianUncorrected
        (
            const surfaceScalarField& gammaMagSf,
            const GeometricField<Type, fvPatchField, volMesh>& vf
        );

        //- Add jump flux to source of matrix
        void addJumpFlux
        (
            fvMatrix<Type>& fvm,
            const surfaceScalarField& gammaMagSf,
            const GeometricField<Type, fvPatchField, volMesh>& vf
        );

// TODO: Two-sided snGradSchemes necessary!
        //- Add correction for both side of face to source of matrix
        void addSnGradsCorrection
        (
            fvMatrix<Type>& fvm,
            const surfaceScalarField& gammaMagSf,
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
                    surfaceInterpolationScheme<GType>::New(mesh, is)
                ),
// TODO: Two-sided snGradSchemes necessary!
                tmp<snGradScheme<Type> >
                (
                    new uncorrectedSnGrad<Type>(mesh)
                )
            ),
// TODO: Two-sided snGradSchemes necessary!
            snGradSchemeName_(is),
            corrected_(true ? snGradSchemeName_ == "corrected" : false),
            limiter_(1.0),
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
            jumpFluxPtr_(NULL),
            vfInterpolScheme_(mesh, gamma_, jumpFluxPtr_)
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

            checkSchemes();
        }


        //- Construct from mesh, gamma and jump flux pointer
        jumpGaussLaplacian
        (
            const fvMesh& mesh,
            const volScalarField& gamma,
            const GeometricField<Type, fvsPatchField, surfaceMesh>* jumpFluxPtr = NULL,
            const word& snGradSchemeName = "corrected"
        )
        :
            laplacianScheme<Type, GType>
            (
                mesh,
                tmp<surfaceInterpolationScheme<GType> >
                (
                    new harmonic<GType>(mesh)
                ),
// TODO: Two-sided snGradSchemes necessary!
                tmp<snGradScheme<Type> >
                (
                    new uncorrectedSnGrad<Type>(mesh)
                )
            ),
// TODO: Two-sided snGradSchemes necessary!
            snGradSchemeName_(snGradSchemeName),
            corrected_(true ? snGradSchemeName_ == "corrected" : false),
            limiter_(1.0),
            gamma_(gamma),
            jumpFluxPtr_(jumpFluxPtr),
            vfInterpolScheme_(mesh, gamma_, jumpFluxPtr_)
        {
            checkSchemes();
        }

    // Destructor

        virtual ~jumpGaussLaplacian()
        {}

    // Member Functions

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
