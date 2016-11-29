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
#include "fvcJumpSurfaceIntegrate.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcSurfaceIntegrate.H"
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
void jumpGaussLaplacian<Type, GType>::checkSchemes() const
{
    word schemeSigma = this->mesh().schemesDict().interpolationScheme
    (
        "interpolate(" + gamma_.name() + ")"
    );

    if
    (
        this->tinterpGammaScheme_().type() != "harmonic"
     || schemeSigma != "harmonic"
    )
    {
        FatalErrorIn
        (
            "jumpGaussLaplacian::jumpGaussLaplacian(...)"
        )   << "For this Laplacian, gamma needs to be"
            << " interpolated harmonically!"
            << exit(FatalError);
    }
}


template<class Type, class GType>
tmp<fvMatrix<Type> >
jumpGaussLaplacian<Type, GType>::fvmLaplacianUncorrected
(
    const surfaceScalarField& gammaMagSf,
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
            gammaMagSf.dimensions()*vf.dimensions()/dimLength
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    fvm.upper() = gammaMagSf.internalField()*deltaCoeffs.internalField();
    fvm.negSumDiag();

    forAll(fvm.psi().boundaryField(), patchI)
    {
        const fvPatchField<Type>& psf = fvm.psi().boundaryField()[patchI];
        const fvsPatchScalarField& patchGamma =
            gammaMagSf.boundaryField()[patchI];

        fvm.internalCoeffs()[patchI] = patchGamma*psf.gradientInternalCoeffs();
        fvm.boundaryCoeffs()[patchI] = -patchGamma*psf.gradientBoundaryCoeffs();
    }

    return tfvm;
}


template<class Type, class GType>
void
jumpGaussLaplacian<Type, GType>::addJumpFlux
(
    fvMatrix<Type>& fvm,
    const surfaceScalarField& gammaMagSf,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    // NOTE: gammaMagSf contains harmonically interpolated gamma!
    const scalarField& gammaMagSfIn = gammaMagSf.internalField();

    const fvMesh& mesh = this->mesh();

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
            dimensioned<Type>
            (
                word(),
                gammaMagSf.dimensions()*vf.dimensions()/dimLength,
                pTraits<Type>::zero
            )
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
            dimensioned<Type>
            (
                word(),
                gammaMagSf.dimensions()*vf.dimensions()/dimLength,
                pTraits<Type>::zero
            )
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& jumpNeiFlux =
        tjumpNeiFlux();
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

    // Secant factors for one-sided gradient
    const surfaceScalarField secAlpha = vfInterpolScheme_.secAlphaOneSided()();
    const scalarField& secAlphaIn = secAlpha.internalField();

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
        // NOTE: The 'mag' is currently only a workaround to avoid
        //       template specializations for different types
        scalar gammaOwn = mag(gammaIn[own]);
        scalar gammaNei = mag(gammaIn[nei]);

        jumpOwnFluxIn[faceI] = secAlphaIn[faceI]
                             * gammaMagSfIn[faceI]/magSfIn[faceI]
                             * (1.0 - gammaOwn/gammaNei) * wP
                             * jumpFluxIn[faceI];

        source[own] -= jumpOwnFluxIn[faceI];

        jumpNeiFluxIn[faceI] = secAlphaIn[faceI]
                             * gammaMagSfIn[faceI]/magSfIn[faceI]
                             * (1.0 - gammaNei/gammaOwn) * wN
                             * jumpFluxIn[faceI];

        source[nei] += jumpNeiFluxIn[faceI];
    }

    // Boundary contributions from jump flux
    // (Contribution will be zero if grad(gamma)*n == 0)
    forAll (mesh.boundary(), patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];

        const scalarField& gammaMagSfPatch = gammaMagSf.boundaryField()[patchI];
        Field<Type>& jumpOwnFluxPatch = jumpOwnFlux.boundaryField()[patchI];

        const unallocLabelList& faceCells = patch.patch().faceCells();

        const scalarField& weightsPatch = weights.boundaryField()[patchI];
        const scalarField& magSfPatch = magSf.boundaryField()[patchI];

        const scalarField& secAlphaPatch = secAlpha.boundaryField()[patchI];

        const Field<GType>& gammaPatch = gamma.boundaryField()[patchI];
        const Field<Type>& jumpFluxPatch = jumpFlux.boundaryField()[patchI];

        forAll (patch, faceI)
        {
            // Cell label
            const label own = faceCells[faceI];

            // Interpolation weights
            scalar wP = weightsPatch[faceI];

            // Cell sigma
            // NOTE: The 'mag' is currently only a workaround to avoid
            //       template specializations for different types
            scalar gammaOwn = mag(gammaIn[own]);
            scalar gammaNei = mag(gammaPatch[faceI]);

            jumpOwnFluxPatch[faceI] = secAlphaPatch[faceI]
                                    * gammaMagSfPatch[faceI]/magSfPatch[faceI]
                                    * (1.0 - gammaOwn/gammaNei) * wP
                                    * jumpFluxPatch[faceI];

            source[own] -= jumpOwnFluxPatch[faceI];
        }
    }

//     // Store jump face fluxes if required
//     if (mesh.schemesDict().fluxRequired(vf.name()))
//     {
//         fvm.faceFluxCorrectionPtr() = tjumpOwnFlux.ptr();
//         fvm.jumpFaceFluxCorrectionPtr() = tjumpNeiFlux.ptr();
//     }
}


template<class Type, class GType>
void jumpGaussLaplacian<Type, GType>::addSnGradsCorrection
(
    fvMatrix<Type>& fvm,
    const surfaceScalarField& gammaMagSf,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    // NOTE: gammaMagSf contains harmonically interpolated gamma!

    const fvMesh& mesh = this->mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tfaceOwnFluxCorr
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "faceOwnFluxCorr("+vf.name()+')',
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>
            (
                word(),
                gammaMagSf.dimensions()*vf.dimensions()/dimLength,
                pTraits<Type>::zero
            )
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& faceOwnFluxCorr =
        tfaceOwnFluxCorr();
    Field<Type>& faceOwnFluxCorrIn = faceOwnFluxCorr.internalField();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tfaceNeiFluxCorr
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "faceNeiFluxCorr("+vf.name()+')',
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>
            (
                word(),
                gammaMagSf.dimensions()*vf.dimensions()/dimLength,
                pTraits<Type>::zero
            )
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& faceNeiFluxCorr =
        tfaceNeiFluxCorr();
    Field<Type>& faceNeiFluxCorrIn = faceNeiFluxCorr.internalField();

    if (!mesh.orthogonal())
    {
Info << "DEBUG: Apply non-orthogonal correction!" << endl;

// TODO: What if jumpFlux not present? We can use normal correction instead...

// TODO: This is hard-coded corrected. How can we use existing snGradSchemes and limiters?

// TODO: What about limiting?

        const Field<Type>& vfIn = vf.internalField();

        // Flux-conserving interpolation of vf
        const GeometricField<Type, fvsPatchField, surfaceMesh> sf =
            vfInterpolScheme_.interpolate(vf);
        const Field<Type>& sfIn = sf.internalField();

        // Cell-based gradient of vf
        typedef typename outerProduct<vector, Type>::type GradType;
        tmp<GeometricField<GradType, fvPatchField, volMesh> > tvfGrad
        (
            fvc::grad(sf)
        );
        tvfGrad().rename("grad(" + vf.name() + ")");

        // Face-interpolated gradient
        tmp<GeometricField<GradType, fvsPatchField, surfaceMesh> > tvfsGrad
        (
            fvc::interpolate(tvfGrad())
        );
        const GeometricField<GradType, fvsPatchField, surfaceMesh>& vfsGrad =
            tvfsGrad();

        tvfGrad.clear();

        // Secant factors for one-sided gradient
        const surfaceScalarField cosAlpha
        (
            1.0/vfInterpolScheme_.secAlphaOneSided()
        );
        const scalarField& cosAlphaIn = cosAlpha.internalField();

        // Mesh addressing
        const unallocLabelList& owner = mesh.owner();
        const unallocLabelList& neighbour = mesh.neighbour();

        // Mesh and basic surface interpolation data
        const surfaceScalarField& weights = mesh.weights();
        const scalarField& weightsIn = weights.internalField();
        const surfaceScalarField& deltaCoeffs = mesh.deltaCoeffs();
        const scalarField& deltaCoeffsIn = deltaCoeffs.internalField();
        const surfaceVectorField& Sf = mesh.Sf();
        const surfaceScalarField& magSf = mesh.magSf();
        const surfaceVectorField& Kf = mesh.correctionVectors();

        // Face-normal vectors
        const surfaceVectorField Nf = Sf/magSf;
        const vectorField NfIn = Nf.internalField();

        // Face-normal part of non-orthogonal correction vector
        const surfaceVectorField KfNf = Nf*(Nf & Kf);
        const vectorField KfNfIn = KfNf.internalField();

        // Face-tangential part of the gradient shows no jump discontinuity
        // and may be corrected from interpolated gradient
        faceOwnFluxCorr = gammaMagSf * (Kf - KfNf) & vfsGrad;
        faceNeiFluxCorr = faceOwnFluxCorr;

        tvfsGrad.clear();

        // Face-normal part of the gradient shows a jump discontinuity such
        // that each side needs to be corrected based on one-sided gradients
        forAll(owner, faceI)
        {
            // Cell labels
            label own = owner[faceI];
            label nei = neighbour[faceI];

            Type snGradOwn = cosAlphaIn[faceI]
                           * (sfIn[faceI] - vfIn[own])
                           / (1.0 - weightsIn[faceI])
                           * deltaCoeffsIn[faceI];

            Type snGradNei = cosAlphaIn[faceI]
                           * (vfIn[nei] - sfIn[faceI])
                           / weightsIn[faceI]
                           * deltaCoeffsIn[faceI];

            faceOwnFluxCorrIn[own] += gammaMagSf[faceI]
                                    * (KfNfIn[faceI] & NfIn[faceI])
                                    * snGradOwn;

            faceNeiFluxCorrIn[nei] += gammaMagSf[faceI]
                                    * (KfNfIn[faceI] & NfIn[faceI])
                                    * snGradNei;
        }

//         // Calculate boundary
//         forAll (mesh.boundary(), patchI)
//         {
//             const fvPatch& patch = mesh.boundary()[patchI];
//
//             const scalarField& gammaMagSfPatch =
//                 gammaMagSf.boundaryField()[patchI];
//             const Field<Type>& sfPatch = sf.boundaryField()[patchI];
//
//             const unallocLabelList& faceCells = patch.patch().faceCells();
//
// //             const vectorField& CfPatch = Cf.boundaryField()[patchI];
//             const scalarField& weightsPatch = weights.boundaryField()[patchI];
//             const scalarField& deltaCoeffsPatch =
//                 deltaCoeffs.boundaryField()[patchI];
//             const vectorField& NfPatch = Nf.boundaryField()[patchI];
//
//             const vectorField& KfNfPatch = KfNf.boundaryField()[patchI];
//
//             const scalarField& cosAlphaPatch = cosAlpha.boundaryField()[patchI];
//
//             Field<Type>& faceOwnFluxCorrPatch =
//                 faceOwnFluxCorr.boundaryField()[patchI];
//
//             forAll (patch, faceI)
//             {
//                 const label own = faceCells[faceI];
//
// //                 Type snGradOwn = cosAlphaPatch[faceI]
// //                                * (sfPatch[faceI] - vfIn[own])
// //                                / mag(CfPatch[faceI] - Cin[own]);
//
//                 Type snGradOwn = cosAlphaPatch[faceI]
//                                * (sfPatch[faceI] - vfIn[own])
//                                / (1.0 - weightsPatch[faceI])
//                                * deltaCoeffsPatch[faceI];
//
//                 faceOwnFluxCorrPatch[own] += gammaMagSfPatch[faceI]
//                                            * (KfNfPatch[faceI] & NfPatch[faceI])
//                                            * snGradOwn;
//             }
//         }

        // Add correction to source
        fvm.source() -=
            mesh.V()*fvc::jumpSurfaceIntegrate<Type>
            (
                faceOwnFluxCorr,
                faceNeiFluxCorr
            )().internalField();

//         // Store face flux corrections if required
//         if (mesh.schemesDict().fluxRequired(vf.name()))
//         {
//             if (!fvm.faceFluxCorrectionPtr())
//             {
//                 fvm.faceFluxCorrectionPtr() = tfaceOwnFluxCorr.ptr();
//             }
//             else
//             {
//                 *fvm.faceFluxCorrectionPtr() += faceOwnFluxCorr;
//                 tfaceOwnFluxCorr.clear();
//             }
//
//             if (!fvm.jumpFaceFluxCorrectionPtr())
//             {
//                 fvm.jumpFaceFluxCorrectionPtr() = tfaceNeiFluxCorr.ptr();
//             }
//             else
//             {
//                 *fvm.jumpFaceFluxCorrectionPtr() += faceNeiFluxCorr;
//                 tfaceNeiFluxCorr.clear();
//             }
//         }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<fvMatrix<Type> >
jumpGaussLaplacian<Type, GType>::fvmLaplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gammaf,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    // NOTE: gammaMagSf contains harmonically interpolated gamma!

    const fvMesh& mesh = this->mesh();

    GeometricField<scalar, fvsPatchField, surfaceMesh> gammaMagSf
    (
        gammaf*mesh.magSf()
    );

    tmp<fvMatrix<Type> > tfvm = fvmLaplacianUncorrected(gammaMagSf, vf);
    fvMatrix<Type>& fvm = tfvm();

    if (jumpFluxPtr_)
    {
        addJumpFlux(fvm, gammaMagSf, vf);
    }

    if (corrected_)
    {
        addSnGradsCorrection(fvm, gammaMagSf, vf);
    }

    return tfvm;
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
    notImplemented(type() + "::fvcLaplacian(vf)");

    const fvMesh& mesh = this->mesh();

    tmp<GeometricField<Type, fvPatchField, volMesh> > tLaplacian
    (
// TODO: This is totally wrong! snGrad does not exist!
        fvc::div(this->tsnGradScheme_().snGrad(vf)*mesh.magSf())
    );

    tLaplacian().rename("laplacian(" + vf.name() + ')');

    return tLaplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
jumpGaussLaplacian<Type, GType>::fvcLaplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gammaf,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    // NOTE: gammaMagSf contains harmonically interpolated gamma!

    notImplemented(type() + "::fvcLaplacian(gammaf, vf)");

    const fvMesh& mesh = this->mesh();

    tmp<GeometricField<Type, fvPatchField, volMesh> > tLaplacian
    (
// TODO: This is totally wrong! snGrad does not exist!
        fvc::div(gammaf*this->tsnGradScheme_().snGrad(vf)*mesh.magSf())
    );

    tLaplacian().rename("laplacian(" + gammaf.name() + ',' + vf.name() + ')');

    return tLaplacian;
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
        vf
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
