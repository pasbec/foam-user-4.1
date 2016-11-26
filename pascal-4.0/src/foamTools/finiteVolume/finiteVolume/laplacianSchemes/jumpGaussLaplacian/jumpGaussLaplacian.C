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

#include "jumpGaussLaplacian.H"
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
jumpGaussLaplacian<Type, GType>::fvmLaplacianUncorrected
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

    // Add jump flux contributions
    if (jumpFluxPtr_)
    {
        fvmLaplacianJumpFlux(fvm, gammaMagSf, vf);
    }

    return tfvm;
}


template<class Type, class GType>
void
jumpGaussLaplacian<Type, GType>::fvmLaplacianJumpFlux
(
    fvMatrix<Type>& fvm,
    const surfaceScalarField& gammaMagSf,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = this->mesh();

    // NOTE: gammaMagSf should contain harmonically interpolated gamma!
    const scalarField& gammaMagSfIn = gammaMagSf.internalField();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tjumpOwnFlux
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "jumpOwnFlux("+vf.name()+')',
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            vf.dimensions()*mesh.deltaCoeffs().dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& jumpOwnFlux =
        tjumpOwnFlux();
    Field<Type>& jumpOwnFluxIn = jumpOwnFlux.internalField();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tjumpNeiFlux
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "jumpNeiFlux("+vf.name()+')',
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            vf.dimensions()*mesh.deltaCoeffs().dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& jumpNeiFlux =
        tjumpOwnFlux();
    Field<Type>& jumpNeiFluxIn = jumpNeiFlux.internalField();

    // Reference to matrix source
    Field<Type>& source = fvm.source();
    // Mesh addressing
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    // Mesh and basic surface interpolation data
    const surfaceScalarField& weights = mesh.weights();
    const scalarField& weightsIn = weights.internalField();
    const surfaceScalarField& magSf = mesh.magSf();
    const scalarField& magSfIn = magSf.internalField();

    // Gamma
    const GeometricField<GType, fvPatchField, volMesh>& gamma = gamma_;
    const Field<GType>& gammaIn = gamma.internalField();

    // Jump flux
    const GeometricField<Type, fvsPatchField, surfaceMesh>& jumpFlux =
        *jumpFluxPtr_;
    const Field<Type>& jumpFluxIn = jumpFlux.internalField();

    // Source contributions from jump flux
    forAll(owner, faceI)
    {
        // Cell labels
        label own = owner[faceI];
        label nei = neighbour[faceI];

        // Interpolation weights
        scalar wP = weightsIn[faceI];
        scalar wN = 1.0 - wP;

        // Cell sigma
// TODO: The 'mag' is currently only a workaround to avoid
//       template specializations for different types
        scalar sigmaOwn = mag(gammaIn[own]);
        scalar sigmaNei = mag(gammaIn[nei]);

        jumpOwnFluxIn[faceI] = gammaMagSfIn[faceI]
                             * (1.0 - sigmaOwn/sigmaNei) * wP
                             * jumpFluxIn[faceI]/magSfIn[faceI];

        source[own] -= jumpOwnFluxIn[faceI];

        jumpNeiFluxIn[faceI] = gammaMagSfIn[faceI]
                             * (1.0 - sigmaNei/sigmaOwn) * wN
                             * jumpFluxIn[faceI]/magSfIn[faceI];

        source[nei] += jumpNeiFluxIn[faceI];
    }

    // Boundary contributions from jump flux
    forAll (mesh.boundary(), patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];

        const scalarField& gammaMagSfPatch = gammaMagSf.boundaryField()[patchI];
        Field<Type>& jumpOwnFluxPatch = jumpOwnFlux.boundaryField()[patchI];

        const unallocLabelList& faceCells = patch.patch().faceCells();

        const scalarField& weightsPatch = weights.boundaryField()[patchI];
        const scalarField& magSfPatch = magSf.boundaryField()[patchI];

        const Field<GType>& gammaPatch = gamma.boundaryField()[patchI];
        const Field<Type>& jumpFluxPatch = jumpFlux.boundaryField()[patchI];

        forAll (patch, faceI)
        {
            // Cell label
            const label own = faceCells[faceI];

            // Interpolation weights
            scalar wP = weightsPatch[faceI];

            // Cell sigma
// TODO: The 'mag' is currently only a workaround to avoid
//       template specializations for different types
            scalar sigmaOwn = mag(gammaIn[own]);
            scalar sigmaNei = mag(gammaPatch[faceI]);

            jumpOwnFluxPatch[faceI] = gammaMagSfPatch[faceI]
                                    * (1.0 - sigmaOwn/sigmaNei) * wP
                                    * jumpFluxPatch[faceI]/magSfPatch[faceI];

            source[own] -= jumpOwnFluxPatch[faceI];
        }
    }

    // Store jump face fluxes if required
    if (mesh.schemesDict().fluxRequired(vf.name()))
    {
        fvm.faceFluxCorrectionPtr() = tjumpOwnFlux.ptr();
        fvm.jumpFaceFluxCorrectionPtr() = tjumpNeiFlux.ptr();
    }
}


template<class Type, class GType>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
jumpGaussLaplacian<Type, GType>::gammaSnGradCorr
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
jumpGaussLaplacian<Type, GType>::fvmLaplacian
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
jumpGaussLaplacian<Type, GType>::fvmLaplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fvmLaplacian(gamma, gamma.name(), vf);
}


template<class Type, class GType>
tmp<fvMatrix<Type> >
jumpGaussLaplacian<Type, GType>::fvmLaplacian
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
jumpGaussLaplacian<Type, GType>::fvcLaplacian
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    FatalErrorIn
    (
        "jumpGaussLaplacian::fvcLaplacian(...)"
    )   << "For this Laplacian, calculus without gamma is disabled."
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
jumpGaussLaplacian<Type, GType>::fvcLaplacian
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
jumpGaussLaplacian<Type, GType>::fvcLaplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fvcLaplacian(gamma, gamma.name(), vf);
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
jumpGaussLaplacian<Type, GType>::fvcLaplacian
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
