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

Description
    Specialisation of electricPotentialGrad for scalars.

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
    const scalarField& w = mesh.weights();
    const scalarField& d = mesh.deltaCoeffs();
    const scalarField& magSf = mesh.magSf();
    const vectorField& Sf = mesh.Sf();

    if (baseScheme() == GAUSS)
    {
        // Internal contributions
        forAll(own, faceI)
        {
            // Cell labels
            label ownFaceI = own[faceI];
            label neiFaceI = nei[faceI];

            // Interpolation weights
            scalar wP = w[faceI];
            scalar wN = 1.0 - wP;

            // Cell sigma
            scalar sigmaOwn = sigmaIn[ownFaceI];
            scalar sigmaNei = sigmaIn[neiFaceI];
            scalar sigmaOff = sigmaNei - sigmaOwn;

            // Face value based on jump conditions
            scalar ssv = sigmaOwn * wP * vsf[ownFaceI]
                       + sigmaNei * wN * vsf[neiFaceI]
// TODO FIXME: deltaCoeffs should come from snGradscheme! Otherwise this is not consistent with
//             the corresponding laplacian operator. Thus, currently only an uncorrected version
//             of the electricPotentialLaplacian exists
                       + sigmaOff * wP * wN * ddtAfluxIn[faceI]/magSf[faceI] / d[faceI];

            ssv /= sigmaMeanIn[faceI];

            // Gradient contributions
            gradIn[ownFaceI] += Sf[faceI] * ssv;
            gradIn[neiFaceI] -= Sf[faceI] * ssv;
        }

        tsigmaMean.clear();

        // Boundary contributions
        // TODO: How can we do this without another interpolation?
        // TODO: Is this enough? What if there are symmetry,
        //       empty or other special bc?
        {
            tmp<surfaceScalarField> tssf
            (
                linearInterpolate(vsf)
            );
            surfaceScalarField& ssf = tssf();

            forAll(mesh.boundary(), patchi)
            {
                const unallocLabelList& pFaceCells =
                    mesh.boundary()[patchi].faceCells();

                const vectorField& pSf = mesh.Sf().boundaryField()[patchi];

                const fvsPatchScalarField& pssf = ssf.boundaryField()[patchi];

                forAll(mesh.boundary()[patchi], facei)
                {
                    grad[pFaceCells[facei]] += pSf[facei]*pssf[facei];
                }
            }

            tssf.clear();
        }

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
