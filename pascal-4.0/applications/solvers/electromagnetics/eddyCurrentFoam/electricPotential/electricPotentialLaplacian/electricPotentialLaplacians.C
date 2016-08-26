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

\*---------------------------------------------------------------------------*/

#include "electricPotentialLaplacian.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    makeFvLaplacianScheme(electricPotentialLaplacian)
}
}


#define declareElectricPotentialLaplacianScalarGamma(Type)                   \
                                                                             \
template<>                                                                   \
Foam::tmp<Foam::fvMatrix<Foam::Type> >                                       \
Foam::fv::electricPotentialLaplacian<Foam::Type, Foam::scalar>::fvmLaplacian \
(                                                                            \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& gamma,         \
    const word& vfGammaName,                                                 \
    const GeometricField<Type, fvPatchField, volMesh>& vf                    \
)                                                                            \
{                                                                            \
    const fvMesh& mesh = this->mesh();                                       \
                                                                             \
    GeometricField<scalar, fvsPatchField, surfaceMesh> gammaMagSf            \
    (                                                                        \
        gamma*mesh.magSf()                                                   \
    );                                                                       \
                                                                             \
    tmp<fvMatrix<Type> > tfvm =                                              \
        fvmLaplacianUncorrected(gammaMagSf, vfGammaName, vf);                \
    fvMatrix<Type>& fvm = tfvm();                                            \
                                                                             \
    if (this->tsnGradScheme_().corrected())                                  \
    {                                                                        \
        if (mesh.schemesDict().fluxRequired(vf.name()))                      \
        {                                                                    \
            fvm.faceFluxCorrectionPtr() =                                    \
                new                                                          \
                GeometricField<Type, fvsPatchField, surfaceMesh>             \
                (                                                            \
                    gammaMagSf*this->tsnGradScheme_().correction(vf)         \
                );                                                           \
                                                                             \
            fvm.source() -=                                                  \
                mesh.V()*                                                    \
                fvc::div                                                     \
                (                                                            \
                    *fvm.faceFluxCorrectionPtr()                             \
                )().internalField();                                         \
        }                                                                    \
        else                                                                 \
        {                                                                    \
            fvm.source() -=                                                  \
                mesh.V()*                                                    \
                fvc::div                                                     \
                (                                                            \
                    gammaMagSf*this->tsnGradScheme_().correction(vf)         \
                )().internalField();                                         \
        }                                                                    \
    }                                                                        \
                                                                             \
    return tfvm;                                                             \
}                                                                            \
                                                                             \
                                                                             \
template<>                                                                   \
Foam::tmp<Foam::fvMatrix<Foam::Type> >                                       \
Foam::fv::electricPotentialLaplacian<Foam::Type, Foam::scalar>::fvmLaplacian \
(                                                                            \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& gamma,         \
    const GeometricField<Type, fvPatchField, volMesh>& vf                    \
)                                                                            \
{                                                                            \
    FatalErrorIn                                                             \
    (                                                                        \
        "electricPotentialLaplacian::fvmLaplacian(...)"                      \
    )   << "Laplacian method with gamma as surfaceField is not allowed."     \
        << exit(FatalError);                                                 \
                                                                             \
    return fvmLaplacian(gamma, word(), vf);                                  \
}                                                                            \
                                                                             \
                                                                             \
template<>                                                                   \
Foam::tmp<Foam::GeometricField<Foam::Type, Foam::fvPatchField, Foam::volMesh> >\
Foam::fv::electricPotentialLaplacian<Foam::Type, Foam::scalar>::fvcLaplacian \
(                                                                            \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& gamma,         \
    const word& vfGammaName,                                                 \
    const GeometricField<Type, fvPatchField, volMesh>& vf                    \
)                                                                            \
{                                                                            \
    const fvMesh& mesh = this->mesh();                                       \
                                                                             \
    tmp<GeometricField<Type, fvPatchField, volMesh> > tLaplacian             \
    (                                                                        \
        fvc::div(gamma*this->tsnGradScheme_().snGrad(vf)*mesh.magSf())       \
    );                                                                       \
                                                                             \
    tLaplacian().rename("laplacian(" + gamma.name() + ',' + vf.name() + ')');\
                                                                             \
    return tLaplacian;                                                       \
}                                                                            \
                                                                             \
                                                                             \
template<>                                                                   \
Foam::tmp<Foam::GeometricField<Foam::Type, Foam::fvPatchField, Foam::volMesh> >\
Foam::fv::electricPotentialLaplacian<Foam::Type, Foam::scalar>::fvcLaplacian \
(                                                                            \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& gamma,         \
    const GeometricField<Type, fvPatchField, volMesh>& vf                    \
)                                                                            \
{                                                                            \
    FatalErrorIn                                                             \
    (                                                                        \
        "electricPotentialLaplacian::fvmLaplacian(...)"                      \
    )   << "Laplacian calculus with gamma as surfaceField is not allowed."   \
        << exit(FatalError);                                                 \
                                                                             \
    return fvcLaplacian(gamma, word(), vf);                                  \
}


declareElectricPotentialLaplacianScalarGamma(scalar);
declareElectricPotentialLaplacianScalarGamma(vector);
declareElectricPotentialLaplacianScalarGamma(sphericalTensor);
declareElectricPotentialLaplacianScalarGamma(symmTensor);
declareElectricPotentialLaplacianScalarGamma(tensor);


// ************************************************************************* //
