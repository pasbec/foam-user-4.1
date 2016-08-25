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

#include "electricPotentialLaplacian.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<fvMatrix<Type> >
electricPotentialLaplacian<Type, GType>::fvmLaplacianUncorrected
(
    const surfaceScalarField& gammaMagSf,
    const word& vfGammaName,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    // NOTE: gammaMagSf contains harmonically interpolated gamma!

    tmp<surfaceScalarField> tdeltaCoeffs =
        this->tsnGradScheme_().deltaCoeffs(vf);
    const surfaceScalarField& deltaCoeffs = tdeltaCoeffs();

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            deltaCoeffs.dimensions()*gammaMagSf.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    fvm.upper() = deltaCoeffs.internalField()*gammaMagSf.internalField();
    fvm.negSumDiag();

    forAll(fvm.psi().boundaryField(), patchI)
    {
        const fvPatchField<Type>& psf = fvm.psi().boundaryField()[patchI];
        const fvsPatchScalarField& patchGamma =
            gammaMagSf.boundaryField()[patchI];

        fvm.internalCoeffs()[patchI] = patchGamma*psf.gradientInternalCoeffs();
        fvm.boundaryCoeffs()[patchI] = -patchGamma*psf.gradientBoundaryCoeffs();
    }

    // Source contributions
    fvm.source() = fvmLaplacianSource(gammaMagSf, vfGammaName);

    return tfvm;
}


template<class Type, class GType>
tmp<Field<Type> >
electricPotentialLaplacian<Type, GType>::fvmLaplacianSource
(
    const surfaceScalarField& gammaMagSf,
    const word& vfGammaName
)
{
    // NOTE: gammaMagSf contains harmonically interpolated gamma!
    const scalarField& gammaMagSfIn = gammaMagSf.internalField();

    const fvMesh& mesh = this->mesh();

    // Electrical conductivity as vol field
    const GeometricField<GType, fvPatchField, volMesh>& sigma =
        mesh.lookupObject<GeometricField<GType, fvPatchField, volMesh> >(vfGammaName);
    const Field<GType>& sigmaIn = sigma.internalField();

    // Rate of change of flux of magnetic vector potential
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ddtAflux =
        mesh.lookupObject<GeometricField<Type, fvsPatchField, surfaceMesh> >(nameDdtAflux_);
    const Field<Type>& ddtAfluxIn = ddtAflux.internalField();

    // Mesh data
    const unallocLabelList& own = mesh.owner();
    const unallocLabelList& nei = mesh.neighbour();
    const surfaceScalarField& w = mesh.weights();
    const scalarField& wIn = w.internalField();
    const surfaceScalarField& magSf = mesh.magSf();
    const scalarField& magSfIn = magSf.internalField();

    tmp<Field<Type> > tsource(new Field<Type>(sigma.size(), pTraits<Type>::zero));
    Field<Type>& source = tsource();


    // Source contributions from jump conditions
    forAll(own, faceI)
    {
        // Cell labels
        label ownFaceI = own[faceI];
        label neiFaceI = nei[faceI];

        // Interpolation weights
        scalar wP = wIn[faceI];
        scalar wN = 1.0 - wP;

        // Cell sigma
// TODO: The 'mag' is currently only a workaround to avoid
//       template specializations for different types
        scalar sigmaOwn = mag(sigmaIn[ownFaceI]);
        scalar sigmaNei = mag(sigmaIn[neiFaceI]);

        source[ownFaceI] -= gammaMagSfIn[faceI]
          * (1.0 - sigmaOwn/sigmaNei) * wP
          * ddtAfluxIn[faceI]/magSfIn[faceI];

        source[neiFaceI] -= gammaMagSfIn[faceI]
          * (sigmaNei/sigmaOwn - 1.0) * wN
          * ddtAfluxIn[faceI]/magSfIn[faceI];
    }

    // Boundary contributions from jump conditions
    forAll (mesh.boundary(), patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];

        // Coupled patches
// TODO: Is this enough? What if there are symmetry, periodic
//       empty or other special bc?
        if (patch.coupled())
        {
            const unallocLabelList& faceCells = patch.patch().faceCells();

            const scalarField& gammaMagSfPatch = gammaMagSf.boundaryField()[patchI];
            const scalarField& wPatch = w.boundaryField()[patchI];
            const scalarField& magSfPatch = magSf.boundaryField()[patchI];
            const Field<GType>& sigmaPatch = sigma.boundaryField()[patchI];
            const Field<Type>& ddtAfluxPatch = ddtAflux.boundaryField()[patchI];

            forAll (patch, faceI)
            {
                // Cell label
                const label ownFaceI = faceCells[faceI];

                // Interpolation weights
                scalar wP = wPatch[faceI];

                // Cell sigma
// TODO: The 'mag' is currently only a workaround to avoid
//       template specializations for different types
                scalar sigmaOwn = mag(sigmaIn[ownFaceI]);
                scalar sigmaNei = mag(sigmaPatch[faceI]);

                source[ownFaceI] -= gammaMagSfPatch[faceI]
                  * (1.0 - sigmaOwn/sigmaNei) * wP
                  * ddtAfluxPatch[faceI]/magSfPatch[faceI];
            }
        }
    }

    return tsource;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
electricPotentialLaplacian<Type, GType>::gammaSnGradCorr
(
    const surfaceVectorField& SfGammaCorr,
    const word& vfGammaName,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = this->mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tgammaSnGradCorr
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "gammaSnGradCorr("+vf.name()+')',
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            SfGammaCorr.dimensions()
           *vf.dimensions()*mesh.deltaCoeffs().dimensions()
        )
    );

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        tgammaSnGradCorr().replace
        (
            cmpt,
            SfGammaCorr & fvc::interpolate(fvc::grad(vf.component(cmpt)))
        );
    }

    return tgammaSnGradCorr;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<fvMatrix<Type> >
electricPotentialLaplacian<Type, GType>::fvmLaplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const word& vfGammaName,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = this->mesh();

    surfaceVectorField Sn = mesh.Sf()/mesh.magSf();

    surfaceVectorField SfGamma = mesh.Sf() & gamma;
    GeometricField<scalar, fvsPatchField, surfaceMesh> SfGammaSn =
        SfGamma & Sn;

    surfaceVectorField SfGammaCorr = SfGamma - SfGammaSn*Sn;

    tmp<fvMatrix<Type> > tfvm =
        fvmLaplacianUncorrected(SfGammaSn, vfGammaName, vf);
    fvMatrix<Type>& fvm = tfvm();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tfaceFluxCorrection
        = gammaSnGradCorr(SfGammaCorr, vfGammaName, vf);

    if (this->tsnGradScheme_().corrected())
    {
        tfaceFluxCorrection() +=
            SfGammaSn*this->tsnGradScheme_().correction(vf);
    }

    fvm.source() -= mesh.V()*fvc::div(tfaceFluxCorrection())().internalField();

    if (mesh.schemesDict().fluxRequired(vf.name()))
    {
        fvm.faceFluxCorrectionPtr() = tfaceFluxCorrection.ptr();
    }

    return tfvm;
}


template<class Type, class GType>
tmp<fvMatrix<Type> >
electricPotentialLaplacian<Type, GType>::fvmLaplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    FatalErrorIn
    (
        "electricPotentialLaplacian::fvmLaplacian(...)"
    )   << "Laplacian method with gamma as surfaceField is not allowed."
        << exit(FatalError);

    // Dummy
    return fvmLaplacian(gamma, word(), vf);
}


template<class Type, class GType>
tmp<fvMatrix<Type> >
electricPotentialLaplacian<Type, GType>::fvmLaplacian
(
    const GeometricField<GType, fvPatchField, volMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fvmLaplacian
    (
        this->tinterpGammaScheme_().interpolate(gamma)(),
        gamma.name(),
        vf
    );
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
electricPotentialLaplacian<Type, GType>::fvcLaplacian
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    FatalErrorIn
    (
        "electricPotentialLaplacian::fvcLaplacian(...)"
    )   << "Laplacian calculus without gamma is not allowed."
        << exit(FatalError);

    const fvMesh& mesh = this->mesh();

    tmp<GeometricField<Type, fvPatchField, volMesh> > tLaplacian
    (
        fvc::div(this->tsnGradScheme_().snGrad(vf)*mesh.magSf())
    );

    tLaplacian().rename("laplacian(" + vf.name() + ')');

    return tLaplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
electricPotentialLaplacian<Type, GType>::fvcLaplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const word& vfGammaName,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = this->mesh();

    surfaceVectorField Sn = mesh.Sf()/mesh.magSf();

    surfaceVectorField SfGamma = mesh.Sf() & gamma;
    GeometricField<scalar, fvsPatchField, surfaceMesh> SfGammaSn
    (
        SfGamma & Sn
    );

    surfaceVectorField SfGammaCorr = SfGamma - SfGammaSn*Sn;

    tmp<GeometricField<Type, fvPatchField, volMesh> > tLaplacian
    (
        fvc::div
        (
            SfGammaSn*this->tsnGradScheme_().snGrad(vf)
          + gammaSnGradCorr(SfGammaCorr, vfGammaName, vf)
        )
    );

    tLaplacian().rename("laplacian(" + gamma.name() + ',' + vf.name() + ')');

    return tLaplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
electricPotentialLaplacian<Type, GType>::fvcLaplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    FatalErrorIn
    (
        "electricPotentialLaplacian::fvmLaplacian(...)"
    )   << "Laplacian calculus with gamma as surfaceField is not allowed."
        << exit(FatalError);

    // Dummy
    return fvcLaplacian(gamma, word(), vf);
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
electricPotentialLaplacian<Type, GType>::fvcLaplacian
(
    const GeometricField<GType, fvPatchField, volMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fvcLaplacian
    (
        this->tinterpGammaScheme_().interpolate(gamma)(),
        gamma.name(),
        vf
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
