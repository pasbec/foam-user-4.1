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
    Basic second-order laplacian using face-gradients and Gauss' theorem.

SourceFiles
    electricPotentialLaplacian.C

\*---------------------------------------------------------------------------*/

#ifndef electricPotentialLaplacian_H
#define electricPotentialLaplacian_H

#include "laplacianScheme.H"
#include "gaussLaplacianScheme.H"

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
    public fv::laplacianScheme<Type, GType>
{

public:

    //- Enumeration defining the valid
    //  base laplacian schemes
    enum baseLaplacianScheme
    {
        GAUSS,
        UNDEFINED
    };

private:

    // Private data

        //- Magnetic properties dictionary
        const IOdictionary& magneticProperties_;

        //- Names of valid base laplacian schemes
        static const NamedEnum<baseLaplacianScheme, 1> baseLaplacianSchemeNames_;

        //- Base laplacian scheme
        baseLaplacianScheme baseLaplacianScheme_;

        //- Name of the faceSet for non-conducting/conducting interface
        word nameInterface_;

        //- Name of the conductivity
        word nameConductivity_;

        //- Name of the electric scalar potential
        word nameElectricPotential_;

        //- Name of the magnetic vector potential
        word nameMagneticPotential_;

    // Private Member Functions

        //- Read dictionary data
        void readMagneticProperties();

        tmp<fvMatrix<Type> > fvmLaplacianUncorrected
        (
            const surfaceScalarField& gammaMagSf,
            GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > gammaSnGradCorr
        (
            const surfaceVectorField& SfGammaCorr,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        //- Disallow default bitwise copy construct
        electricPotentialLaplacian(const electricPotentialLaplacian&);

        //- Disallow default bitwise assignment
        void operator=(const electricPotentialLaplacian&);


public:

    //- Runtime type information
    TypeName("electricPotential");


    // Constructors

        //- Construct null
        electricPotentialLaplacian(const fvMesh& mesh)
        :
            laplacianScheme<Type, GType>(mesh),
            magneticProperties_
            (
                mesh.lookupObject<IOdictionary>
                (
                    "magneticProperties"
                )
            ),
            baseLaplacianScheme_(GAUSS)
        {
            readMagneticProperties();
        }

        //- Construct from Istream
        electricPotentialLaplacian(const fvMesh& mesh, Istream& is)
        :
            laplacianScheme<Type, GType>(mesh, is),
            magneticProperties_
            (
                mesh.lookupObject<IOdictionary>
                (
                    "magneticProperties"
                )
            ),
            baseLaplacianScheme_(UNDEFINED)
        {
            if (is.eof())
            {
                FatalIOErrorIn
                (
                    "electricPotentialLaplacian<Type, GType>::electricPotentialLaplacian()",
                    is
                )   << "Base laplacian scheme not specified!" << endl
                    << "Valid base gradient schemes are :" << endl
                    << baseLaplacianSchemeNames_.sortedToc()
                    << exit(FatalIOError);
            }

            word nameBaseLaplacianScheme = word(is);

            typename NamedEnum<baseLaplacianScheme, 1>::const_iterator
                cIter = baseLaplacianSchemeNames_.find(nameBaseLaplacianScheme);

            if (cIter == baseLaplacianSchemeNames_.end())
            {
                FatalIOErrorIn
                (
                    "electricPotentialLaplacian<Type, GType>::electricPotentialLaplacian()",
                    is
                )   << "Base gradient scheme " << nameBaseLaplacianScheme
                    << " is not supported!" << endl
                    << "Valid base gradient schemes are :" << endl
                    << baseLaplacianSchemeNames_.sortedToc()
                    << exit(FatalIOError);
            }

            baseLaplacianScheme_ = baseLaplacianSchemeNames_[nameBaseLaplacianScheme];

            readMagneticProperties();
        }

        //- Construct from mesh, interpolation and snGrad schemes
        electricPotentialLaplacian
        (
            const fvMesh& mesh,
            const tmp<surfaceInterpolationScheme<GType> >& igs,
            const tmp<snGradScheme<Type> >& sngs
        )
        :
            laplacianScheme<Type, GType>(mesh, igs, sngs),
            magneticProperties_
            (
                mesh.lookupObject<IOdictionary>
                (
                    "magneticProperties"
                )
            ),
            baseLaplacianScheme_(GAUSS)
        {
            readMagneticProperties();
        }

        //- Construct from mesh, interpolation, snGrad and baseLaplcian schemes
        electricPotentialLaplacian
        (
            const fvMesh& mesh,
            const tmp<surfaceInterpolationScheme<GType> >& igs,
            const tmp<snGradScheme<Type> >& sngs,
            const baseLaplacianScheme bls
        )
        :
            laplacianScheme<Type, GType>(mesh, igs, sngs),
            magneticProperties_
            (
                mesh.lookupObject<IOdictionary>
                (
                    "magneticProperties"
                )
            ),
            baseLaplacianScheme_(bls)
        {
            if (bls == UNDEFINED)
            {
                FatalErrorIn
                (
                    "electricPotentialLaplacian<Type, GType>::electricPotentialLaplacian()"
                )   << "Base laplacian scheme not specified!" << endl
                    << "Valid base gradient schemes are :" << endl
                    << baseLaplacianSchemeNames_.sortedToc()
                    << exit(FatalIOError);
            }

            readMagneticProperties();
        }


    // Destructor

        virtual ~electricPotentialLaplacian()
        {}


    // Member Functions

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcLaplacian
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type> > fvmLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );
};


// Use macros to emulate partial-specialisation of the the Laplacian functions
// for scalar diffusivity gamma

#define defineLaplacianScalarGamma(Type)                                    \
                                                                            \
template<>                                                                  \
tmp<fvMatrix<Type> > electricPotentialLaplacian<Type, scalar>::fvmLaplacian \
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
