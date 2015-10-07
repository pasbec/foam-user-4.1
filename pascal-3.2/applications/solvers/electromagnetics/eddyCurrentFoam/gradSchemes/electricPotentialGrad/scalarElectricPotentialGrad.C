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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


template<>
tmp<GeometricField<vector, fvPatchField, volMesh> >
electricPotentialGrad<scalar>::gradf
(
    const GeometricField<scalar, fvPatchField, volMesh>& vsf,
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& ssf,
    const word& name
) const
{
    const fvMesh& mesh = ssf.mesh();

    tmp<volVectorField> tgGrad
    (
        new volVectorField
        (
            IOobject
            (
                name,
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<vector>
            (
                "0",
                ssf.dimensions()/dimLength,
                pTraits<vector>::zero
            ),
            zeroGradientFvPatchField<vector>::typeName
        )
    );
    volVectorField& gGrad = tgGrad();

    const unallocLabelList& own = mesh.owner();
    const unallocLabelList& nei = mesh.neighbour();

    const vectorField& Sf = mesh.Sf();
    const scalarField& magSf = mesh.magSf();
    const scalarField& w = mesh.weights();
    const scalarField& d = mesh.deltaCoeffs();

    // TODO: Optimize speed with tmp
    const vectorField Sn = Sf/magSf;

    vectorField& igGrad = gGrad;
    const scalarField& issf = ssf;

    forAll(own, facei)
    {
        vector Sfssf = Sf[facei]*issf[facei];

        igGrad[own[facei]] += Sfssf;
        igGrad[nei[facei]] -= Sfssf;
    }

    // TODO: Optimize speed with tmp
    const labelList cnFaces =
        mesh.lookupObject<faceSet>
        (
            nameInterface_
        ).toc();

    const dimensionedScalar frequency
    (
        magneticProperties_.lookup("frequency")
    );

    const dimensionedScalar omega
    (
        "omega",
        2.0 * mathematicalConstant::pi * frequency
    );

    // TODO: Prevent segfault!
    word complexNameV = vsf.name();
    word complexPartV =
        complexNameV
        (
            complexNameV.size()-2,
            complexNameV.size()-1
        );

    word complexPartA = word();
    scalar complexSignA = 0.0;

    if (complexPartV == "Re")
    {
        complexPartA = "Im";
        complexSignA = 1.0;
    }
    else if (complexPartV == "Im")
    {
        complexPartA = "Re";
        complexSignA = -1.0;
    }
    else
    {
        // TODO: Error handling!
    }

    const volVectorField& A =
        mesh.lookupObject<volVectorField>
        (
            nameMagneticPotential_ + complexPartA
        );

    scalar As = complexSignA;

    // TODO: Optimize speed with tmp SfsA
    forAll(cnFaces, facei)
    {
        label faceI = cnFaces[facei];

        // Remove wrongly interpolated gradient
        // contributions from cells at the interface
        // between non-conducting and conducting region
        {
            vector Sfssf = Sf[faceI]*issf[faceI];

            igGrad[own[faceI]] -= Sfssf;
            igGrad[nei[faceI]] += Sfssf;
        }

        // Add correct face value contributions to cells
        // for the electric scalar potential from its
        // corresponding gradient at the interface between
        // non-conducting and conducting region depending on
        // the magnetic vector potential
        {
            // Face centre to cell centre distances
            scalar fN = w[faceI] / d[faceI];
            scalar Pf = 1.0 / d[faceI] - fN;

            // Gradient weight coeffs
            scalar wP = omega.value() * w[faceI];
            scalar wN = omega.value() - wP;

            // Weighted face gradients
            scalar SfAwOwn = wP * As * (A[own[faceI]] & Sn[faceI]);
            scalar SfAwNei = wN * As * (A[nei[faceI]] & Sn[faceI]);

            // Extrapolated face values
            scalar ssfOwn = vsf[own[faceI]] + Pf * (SfAwOwn + SfAwNei);
            scalar ssfNei = vsf[nei[faceI]] - fN * (SfAwOwn + SfAwNei);

            igGrad[own[faceI]] += Sf[faceI] * ssfOwn;
            igGrad[nei[faceI]] -= Sf[faceI] * ssfNei;
        }
    }

    // NOTE: Currently it is NOT SUPPORTED or even tested to put the
    //       conducting region on any type of boundary

    forAll(mesh.boundary(), patchi)
    {
        const unallocLabelList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const vectorField& pSf = mesh.Sf().boundaryField()[patchi];

        const fvsPatchField<scalar>& pssf = ssf.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            igGrad[pFaceCells[facei]] += pSf[facei]*pssf[facei];
        }
    }

    igGrad /= mesh.V();

    gGrad.correctBoundaryConditions();

    return tgGrad;
}


template<>
tmp<GeometricField<vector, fvPatchField, volMesh> >
electricPotentialGrad<scalar>::calcGrad
(
    const GeometricField<scalar, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    tmp<volVectorField> tgGrad
    (
        gradf(vsf, tinterpScheme_().interpolate(vsf), name)
    );
    volVectorField& gGrad = tgGrad();

    gGrad.rename("grad(" + vsf.name() + ')');
    this->correctBoundaryConditions(vsf, gGrad);

    return tgGrad;
}


// // TODO
template<>
tmp<BlockLduSystem<vector, vector> >
electricPotentialGrad<scalar>::fvmGrad
(
    const GeometricField<scalar, fvPatchField, volMesh>& vf
) const
{
    FatalErrorIn
    (
        "tmp<BlockLduSystem> fvmGrad\n"
        "(\n"
        "    GeometricField<Type, fvPatchField, volMesh>&\n"
        ")\n"
    )   << "Implicit gradient operator not yet implemented."
        << abort(FatalError);

    tmp<BlockLduSystem<vector, vector> > tbs
    (
        new BlockLduSystem<vector, vector>(vf.mesh())
    );
    return tbs;
}
// template<>
// tmp<BlockLduSystem<vector, vector> >
// electricPotentialGrad<scalar>::fvmGrad
// (
//     const GeometricField<scalar, fvPatchField, volMesh>& vf
// ) const
// {
//     tmp<surfaceScalarField> tweights = this->tinterpScheme_().weights(vf);
//     const scalarField& wIn = tweights().internalField();
//
//     const fvMesh& mesh = vf.mesh();
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
//         const fvsPatchScalarField& pw = tweights().boundaryField()[patchI];
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
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
