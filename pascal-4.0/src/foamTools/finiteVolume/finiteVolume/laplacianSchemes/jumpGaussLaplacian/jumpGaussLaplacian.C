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
    const volTypeField& vf
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
    const volTypeField& vf
)
{
    // NOTE: gammaMagSf contains harmonically interpolated gamma!
    const scalarField& gammaMagSfIn = gammaMagSf.internalField();

    const fvMesh& mesh = this->mesh();

    tmp<surfaceTypeField> tjumpOwnFlux
    (
        new surfaceTypeField
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
    surfaceTypeField& jumpOwnFlux = tjumpOwnFlux();
    Field<Type>& jumpOwnFluxIn = jumpOwnFlux.internalField();

    tmp<surfaceTypeField> tjumpNeiFlux
    (
        new surfaceTypeField
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
    surfaceTypeField& jumpNeiFlux = tjumpNeiFlux();
    Field<Type>& jumpNeiFluxIn = jumpNeiFlux.internalField();

    // Mesh addressing
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    // Mesh and basic surface interpolation data
    const surfaceScalarField& weights = mesh.weights();
    const scalarField& weightsIn = weights.internalField();
    const surfaceScalarField& magSf = mesh.magSf();
    const scalarField& magSfIn = magSf.internalField();

    // Secant factors for one-sided gradient
    const surfaceScalarField secAlpha = vfFluxConsScheme_.secAlphaOneSided()();
    const scalarField& secAlphaIn = secAlpha.internalField();

    // Gamma
    const volGTypeField& gamma = gamma_;
    const Field<GType>& gammaIn = gamma.internalField();

    // Jump flux
    const surfaceTypeField& jumpFlux = *jumpFluxPtr_;
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

        jumpNeiFluxIn[faceI] = secAlphaIn[faceI]
                             * gammaMagSfIn[faceI]/magSfIn[faceI]
                             * (1.0 - gammaNei/gammaOwn) * wN
                             * jumpFluxIn[faceI];
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
        }
    }

    // Add jump face flux to source
    fvm.source() =
        -mesh.V()*fvc::surfaceIntegrate<Type>
        (
            jumpOwnFlux,
            jumpNeiFlux
        )().internalField();

    // Store jump face fluxes if required
    if (mesh.schemesDict().fluxRequired(vf.name()))
    {
        fvm.faceFluxCorrectionPtr() = tjumpOwnFlux.ptr();
        fvm.jumpFaceFluxCorrectionPtr() = tjumpNeiFlux.ptr();
    }
}


template<class Type, class GType>
void jumpGaussLaplacian<Type, GType>::addSnGradsCorrection
(
    fvMatrix<Type>& fvm,
    const surfaceScalarField& gammaMagSf,
    const volTypeField& vf
)
{
    // NOTE: gammaMagSf contains harmonically interpolated gamma!

    const fvMesh& mesh = this->mesh();

    tmp<surfaceTypeField > tfaceOwnFluxCorr
    (
        new surfaceTypeField
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
    surfaceTypeField& faceOwnFluxCorr = tfaceOwnFluxCorr();
    Field<Type>& faceOwnFluxCorrIn = faceOwnFluxCorr.internalField();

    tmp<surfaceTypeField> tfaceNeiFluxCorr
    (
        new surfaceTypeField
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
    surfaceTypeField& faceNeiFluxCorr = tfaceNeiFluxCorr();
    Field<Type>& faceNeiFluxCorrIn = faceNeiFluxCorr.internalField();

    if (!mesh.orthogonal())
    {
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
        const vectorField& NfIn = Nf.internalField();

        // Face-normal part of non-orthogonal correction vector
        const surfaceVectorField KfNf = Nf*(Nf & Kf);
        const vectorField& KfNfIn = KfNf.internalField();

        // Face-tangential part of the gradient shows no jump discontinuity
        // and may be corrected from interpolated gradient
        tmp<surfaceTypeField> tcorr = this->tsnGradScheme_().correction(vf);
        faceOwnFluxCorr =
            gammaMagSf * tcorr()
          * Kf/(mag(Kf) + VSMALL) & (Kf - KfNf)/(mag(Kf - KfNf) + VSMALL);
        faceNeiFluxCorr = faceOwnFluxCorr;

        tcorr.clear();

        const Field<Type>& vfIn = vf.internalField();

        // Flux-conserving interpolation of vf
        const surfaceTypeField sf = vfFluxConsScheme_.interpolate(vf);
        const Field<Type>& sfIn = sf.internalField();

        // Secant factors for one-sided gradient
        const surfaceScalarField cosAlpha
        (
            1.0/vfFluxConsScheme_.secAlphaOneSided()
        );
        const scalarField& cosAlphaIn = cosAlpha.internalField();

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

            faceOwnFluxCorrIn[faceI] += gammaMagSf[faceI]
                                      * (KfNfIn[faceI] & NfIn[faceI])
                                      * snGradOwn;

            faceNeiFluxCorrIn[faceI] += gammaMagSf[faceI]
                                      * (KfNfIn[faceI] & NfIn[faceI])
                                      * snGradNei;
        }

        forAll (mesh.boundary(), patchI)
        {
            const fvPatch& patch = mesh.boundary()[patchI];

            if (patch.coupled())
            {
                const scalarField& gammaMagSfPatch =
                    gammaMagSf.boundaryField()[patchI];
                const Field<Type>& sfPatch = sf.boundaryField()[patchI];

                const unallocLabelList& faceCells = patch.patch().faceCells();

                const scalarField& weightsPatch =
                    weights.boundaryField()[patchI];
                const scalarField& deltaCoeffsPatch =
                    deltaCoeffs.boundaryField()[patchI];
                const vectorField& NfPatch = Nf.boundaryField()[patchI];

                const vectorField& KfNfPatch = KfNf.boundaryField()[patchI];

                const scalarField& cosAlphaPatch =
                    cosAlpha.boundaryField()[patchI];

                Field<Type>& faceOwnFluxCorrPatch =
                    faceOwnFluxCorr.boundaryField()[patchI];

                forAll (patch, faceI)
                {
                    const label own = faceCells[faceI];

                    Type snGradOwn = cosAlphaPatch[faceI]
                                   * (sfPatch[faceI] - vfIn[own])
                                   / (1.0 - weightsPatch[faceI])
                                   * deltaCoeffsPatch[faceI];

                    faceOwnFluxCorrPatch[faceI] += gammaMagSfPatch[faceI]
                                                 * (KfNfPatch[faceI] & NfPatch[faceI])
                                                 * snGradOwn;
                }
            }
        }

        // Add correction to source
        fvm.source() -=
            mesh.V()*fvc::surfaceIntegrate<Type>
            (
                faceOwnFluxCorr,
                faceNeiFluxCorr
            )().internalField();

        // Store face flux corrections if required
        if (mesh.schemesDict().fluxRequired(vf.name()))
        {
            if (!fvm.faceFluxCorrectionPtr())
            {
                fvm.faceFluxCorrectionPtr() = tfaceOwnFluxCorr.ptr();
            }
            else
            {
                *fvm.faceFluxCorrectionPtr() += faceOwnFluxCorr;
                tfaceOwnFluxCorr.clear();
            }

            if (!fvm.jumpFaceFluxCorrectionPtr())
            {
                fvm.jumpFaceFluxCorrectionPtr() = tfaceNeiFluxCorr.ptr();
            }
            else
            {
                *fvm.jumpFaceFluxCorrectionPtr() += faceNeiFluxCorr;
                tfaceNeiFluxCorr.clear();
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<fvMatrix<Type> >
jumpGaussLaplacian<Type, GType>::fvmLaplacian
(
    const surfaceGTypeField& gammaf,
    const volTypeField& vf
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

    if (this->tsnGradScheme_().corrected())
    {
        addSnGradsCorrection(fvm, gammaMagSf, vf);
    }

    return tfvm;
}


template<class Type, class GType>
tmp<fvMatrix<Type> >
jumpGaussLaplacian<Type, GType>::fvmLaplacian
(
    const volGTypeField& gamma,
    const volTypeField& vf
)
{
    return fvmLaplacian
    (
        this->tinterpGammaScheme_().interpolate(gamma)(),
        vf
    );
}


// TODO
template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
jumpGaussLaplacian<Type, GType>::fvcLaplacian
(
    const volTypeField& vf
)
{
    notImplemented(type() + "::fvcLaplacian(vf)");

    const fvMesh& mesh = this->mesh();

    tmp<volTypeField > tLaplacian
    (
// TODO: Two-sided snGradSchemes necessary!
        fvc::div(this->tsnGradScheme_().snGrad(vf)*mesh.magSf())
    );

    tLaplacian().rename("laplacian(" + vf.name() + ')');

    return tLaplacian;
}


// TODO
template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
jumpGaussLaplacian<Type, GType>::fvcLaplacian
(
    const surfaceGTypeField& gammaf,
    const volTypeField& vf
)
{
    // NOTE: gammaMagSf contains harmonically interpolated gamma!

    const fvMesh& mesh = this->mesh();

    tmp<fvMatrix<Type> > tfvm = fvmLaplacian(gammaf, vf);
    fvMatrix<Type>& fvm = tfvm();

    tmp<volTypeField > tLaplacian
    (
        fvm.A()*vf - fvm.H()
    );

//     tmp<volTypeField > tLaplacian
//     (
//         new volTypeField
//         (
//             IOobject
//             (
//                 "laplacian("+gammaf.name()+","+vf.name()+')',
//                 vf.instance(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE
//             ),
//             mesh,
//             dimensioned<Type>
//             (
//                 word(),
//                 gammaf.dimensions()*vf.dimensions()/dimArea,
//                 pTraits<Type>::zero
//             ),
//             zeroGradientFvPatchField<Type>::typeName
//         )
//     );
//
//     tLaplacian().internalField() = fvm.residual()/mesh.V();
//     tfvm.clear();
//
//     tLaplacian().correctBoundaryConditions();

//     tmp<volTypeField > tLaplacian
//     (
// // TODO: Two-sided snGradSchemes necessary!
//         fvc::div(gammaf*this->tsnGradScheme_().snGrad(vf)*mesh.magSf())
//     );

    tLaplacian().rename("laplacian(" + gammaf.name() + ',' + vf.name() + ')');

    return tLaplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
jumpGaussLaplacian<Type, GType>::fvcLaplacian
(
    const volGTypeField& gamma,
    const volTypeField& vf
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
