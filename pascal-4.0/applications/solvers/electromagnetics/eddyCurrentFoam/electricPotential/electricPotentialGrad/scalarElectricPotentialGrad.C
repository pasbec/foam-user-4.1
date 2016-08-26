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

#include "electricPotentialGrad.H"
#include "gaussGrad.H"
#include "leastSquaresGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<GeometricField<vector, fvPatchField, volMesh> >
electricPotentialGrad<scalar>::calcGrad
(
    const GeometricField<scalar, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    const fvMesh& mesh = vsf.mesh();

    const scalarField& vsfIn = vsf.internalField();

    // Prepare gradient field
    tmp<volVectorField> tgrad
    (
        new volVectorField
        (
            IOobject
            (
                name,
                vsf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<vector>
            (
                "zero",
                vsf.dimensions()/dimLength,
                pTraits<vector>::zero
            ),
            zeroGradientFvPatchField<vector>::typeName
        )
    );
    volVectorField& grad = tgrad();
    vectorField& gradIn = grad.internalField();

    // Electrical conductivity from gamma
    const volScalarField& sigma =
        mesh.lookupObject<volScalarField>(nameSigma_);
    const scalarField& sigmaIn = sigma.internalField();

    // Linear interpolated sigma (mean)
    tmp<surfaceScalarField> tsigmaMean
    (
        linearInterpolate(sigma)
    );
    const surfaceScalarField& sigmaMean = tsigmaMean();
    const scalarField& sigmaMeanIn = sigmaMean.internalField();

    // Rate of change of flux of magnetic vector potential
    const surfaceScalarField& ddtAflux =
        mesh.lookupObject<surfaceScalarField>(nameDdtAflux_);
    const scalarField& ddtAfluxIn = ddtAflux.internalField();

    // Mesh data
    const unallocLabelList& own = mesh.owner();
    const unallocLabelList& nei = mesh.neighbour();
    const surfaceScalarField& w = mesh.weights();
    const scalarField& wIn = w.internalField();
    const surfaceScalarField& d = mesh.deltaCoeffs();
    const scalarField& dIn = d.internalField();
    const surfaceScalarField& magSf = mesh.magSf();
    const scalarField& magSfIn = magSf.internalField();
    const surfaceVectorField& Sf = mesh.Sf();
    const vectorField& SfIn = Sf.internalField();

    if (baseScheme() == GAUSS)
    {
        // Internal contributions
        forAll(own, faceI)
        {
            // Cell labels
            label ownFaceI = own[faceI];
            label neiFaceI = nei[faceI];

            // Interpolation weights
            scalar wP = wIn[faceI];
            scalar wN = 1.0 - wP;

            // Cell sigma
            scalar sigmaOwn = sigmaIn[ownFaceI];
            scalar sigmaNei = sigmaIn[neiFaceI];
            scalar sigmaOff = sigmaNei - sigmaOwn;

            // Face value based on jump conditions
            scalar ssv = sigmaOwn * wP * vsfIn[ownFaceI]
                       + sigmaNei * wN * vsfIn[neiFaceI]
// TODO FIXME: deltaCoeffs should come from snGradscheme! Otherwise this is not consistent with
//             the corresponding laplacian operator. Thus, currently only an uncorrected version
//             of the electricPotentialLaplacian exists
                       + sigmaOff * wP * wN * ddtAfluxIn[faceI]/magSfIn[faceI] / dIn[faceI];

            ssv /= sigmaMeanIn[faceI];

            // Gradient contributions
            gradIn[ownFaceI] += SfIn[faceI] * ssv;
            gradIn[neiFaceI] -= SfIn[faceI] * ssv;
        }

        // Boundary contributions from jump conditions
        forAll (mesh.boundary(), patchI)
        {
            const fvPatch& patch = mesh.boundary()[patchI];

            const unallocLabelList& faceCells = patch.patch().faceCells();

            const scalarField& vsfPatch = vsf.boundaryField()[patchI];

            const vectorField& SfPatch = Sf.boundaryField()[patchI];

            // Coupled patches
// TODO: Is this enough? What if there are symmetry, periodic
//       empty or other special bc?
            if (patch.coupled())
            {
                const scalarField& wPatch = w.boundaryField()[patchI];
                const scalarField& sigmaPatch = sigma.boundaryField()[patchI];
                const scalarField& ddtAfluxPatch = ddtAflux.boundaryField()[patchI];
                const scalarField& magSfPatch = magSf.boundaryField()[patchI];
                const scalarField& dPatch = d.boundaryField()[patchI];
                const scalarField& sigmaMeanPatch = sigmaMean.boundaryField()[patchI];

                forAll (patch, faceI)
                {
                    // Cell label
                    const label ownFaceI = faceCells[faceI];

                    // Interpolation weights
                    scalar wP = wPatch[faceI];
                    scalar wN = 1.0 - wP;

                    // Cell sigma
                    scalar sigmaOwn = sigmaIn[ownFaceI];
                    scalar sigmaNei = sigmaPatch[faceI];
                    scalar sigmaOff = sigmaNei - sigmaOwn;

                    // Face value based on jump conditions
                    scalar ssv = sigmaOwn * wP * vsfIn[ownFaceI]
                               + sigmaNei * wN * vsfPatch[faceI]
// TODO FIXME: deltaCoeffs should come from snGradscheme! Otherwise this is not consistent with
//             the corresponding laplacian operator. Thus, currently only an uncorrected version
//             of the electricPotentialLaplacian exists
                               + sigmaOff * wP * wN * ddtAfluxPatch[faceI]/magSfPatch[faceI] / dPatch[faceI];

                    ssv /= sigmaMeanPatch[faceI];

                    // Gradient contributions
                    gradIn[ownFaceI] += SfPatch[faceI] * ssv;
                }
            }
            else
            {
                forAll (patch, faceI)
                {
                    // Cell label
                    const label ownFaceI = faceCells[faceI];

                    // Gradient contributions
                    gradIn[ownFaceI] += SfPatch[faceI] * vsfPatch[faceI];

                }
            }
        }

        tsigmaMean.clear();

        gradIn /= mesh.V();
    }
    else if (baseScheme() == LEASTSQUARES)
    {
        FatalErrorIn
        (
            "electricPotentialGrad<scalar>::calcGrad(...)"
        )   << "Gradient calculation for LEASTSQUARES "
            << "base scheme not yet implemented."
            << abort(FatalError);
    }

    grad.rename("grad(" + vsf.name() + ')');
    grad.correctBoundaryConditions();

    return tgrad;
}


// TODO: Implicit gradient
template<>
tmp<BlockLduSystem<vector, vector> >
electricPotentialGrad<scalar>::fvmGrad
(
    const GeometricField<scalar, fvPatchField, volMesh>& vf
) const
{
    FatalErrorIn
    (
        "electricPotentialGrad<scalar>::fvmGrad(...)"
    )   << "Implicit gradient operator not yet implemented."
        << abort(FatalError);

    tmp<BlockLduSystem<vector, vector> > tbs
    (
        new BlockLduSystem<vector, vector>(vf.mesh())
    );

    return tbs;

//     const fvMesh& mesh = vf.mesh();
//
//     const surfaceScalarField& weights = mesh.weights();
//     const scalarField& wIn = weights.internalField();
//
//     tmp<BlockLduSystem<vector, vector> > tbs
//     (
//         new BlockLduSystem<vector, vector>(mesh)
//     );
//     BlockLduSystem<vector, vector>& bs = tbs();
//     vectorField& source = bs.source();
//
//     // Grab ldu parts of block matrix as linear always
//     CoeffField<vector>::linearTypeField& d = bs.diag().asLinear();
//     CoeffField<vector>::linearTypeField& u = bs.upper().asLinear();
//     CoeffField<vector>::linearTypeField& l = bs.lower().asLinear();
//
//     const vectorField& SfIn = mesh.Sf().internalField();
//
//     l = -wIn*SfIn;
//     u = l + SfIn;
//     bs.negSumDiag();
//
//     // Boundary contributions
//     forAll (vf.boundaryField(), patchI)
//     {
//         const fvPatchScalarField& pf = vf.boundaryField()[patchI];
//         const fvPatch& patch = pf.patch();
//         const vectorField& pSf = patch.Sf();
//         const fvsPatchScalarField& pw = weights.boundaryField()[patchI];
//         const labelList& fc = patch.faceCells();
//
//         const scalarField internalCoeffs(pf.valueInternalCoeffs(pw));
//
//         // Diag contribution
//         forAll (pf, faceI)
//         {
//             d[fc[faceI]] += internalCoeffs[faceI]*pSf[faceI];
//         }
//
//         if (patch.coupled())
//         {
//             CoeffField<vector>::linearTypeField& pcoupleUpper =
//                 bs.coupleUpper()[patchI].asLinear();
//             CoeffField<vector>::linearTypeField& pcoupleLower =
//                 bs.coupleLower()[patchI].asLinear();
//
//             const vectorField pcl = -pw*pSf;
//             const vectorField pcu = pcl + pSf;
//
//             // Coupling  contributions
//             pcoupleLower -= pcl;
//             pcoupleUpper -= pcu;
//         }
//         else
//         {
//             const scalarField boundaryCoeffs(pf.valueBoundaryCoeffs(pw));
//
//             // Boundary contribution
//             forAll (pf, faceI)
//             {
//                 source[fc[faceI]] -= boundaryCoeffs[faceI]*pSf[faceI];
//             }
//         }
//     }
//
//     // Interpolation schemes with corrections not accounted for
//
//     return tbs;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
