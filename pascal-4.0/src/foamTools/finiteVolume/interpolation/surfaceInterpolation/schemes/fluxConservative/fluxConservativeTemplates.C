
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

#include "fluxConservative.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::surfaceScalarField>
Foam::fluxConservative<Type>::secAlphaNonOrth() const
{
    const fvMesh& mesh = this->mesh();

    tmp<surfaceScalarField> tsecAlpha
    (
        new surfaceScalarField
        (
            IOobject
            (
                "secAlphaNonOrthAngle",
                mesh.time().timeName(),
                mesh
            ),
            mesh ,
            dimensionedScalar(word(), dimless, 1.0)
        )
    );
    surfaceScalarField& secAlpha = tsecAlpha();
    scalarField& secAlphaIn = secAlpha.internalField();

    if (!mesh.orthogonal())
    {
        // Mesh addressing
        const unallocLabelList& owner = mesh.owner();

        // Mesh and basic surface interpolation data
        const volVectorField& C = mesh.C();
        const surfaceScalarField& deltaCoeffs = mesh.deltaCoeffs();
        const surfaceVectorField& Sf = mesh.Sf();
        const vectorField& SfIn = Sf.internalField();
        const surfaceScalarField& magSf = mesh.magSf();
        const scalarField& magSfIn = magSf.internalField();
        const surfaceVectorField& Kf = mesh.correctionVectors();
        const surfaceVectorField Df = Sf/magSf - Kf;
        const vectorField& DfIn = Df.internalField();

        // Calculate internal secAlpha
        forAll (owner, faceI)
        {
            vector NfInI = SfIn[faceI]/magSfIn[faceI];

            secAlphaIn[faceI] = 1.0/(NfInI & DfIn[faceI]);
        }

        // Calculate boundary secAlpha
        forAll (mesh.boundary(), patchI)
        {
            const fvPatch& patch = mesh.boundary()[patchI];

            if (patch.coupled())
            {
                const scalarField& magSfPatch = magSf.boundaryField()[patchI];
                const vectorField& SfPatch = Sf.boundaryField()[patchI];
                const vectorField& DfPatch = Df.boundaryField()[patchI];

                scalarField& secAlphaPatch = secAlpha.boundaryField()[patchI];

                forAll (patch, faceI)
                {
                    vector NfPatchI = SfPatch[faceI]/magSfPatch[faceI];

                    secAlphaPatch[faceI] = 1.0/(NfPatchI & DfPatch[faceI]);
                }
            }
            else
            {
                const unallocLabelList& faceCells = patch.patch().faceCells();

                const scalarField& deltaCoeffsPatch =
                    deltaCoeffs.boundaryField()[patchI];
                const scalarField& magSfPatch = magSf.boundaryField()[patchI];
                const vectorField& SfPatch = Sf.boundaryField()[patchI];
                const vectorField& CfPatch = patch.Cf();

                scalarField& secAlphaPatch = secAlpha.boundaryField()[patchI];

                forAll (patch, faceI)
                {
                    const label own = faceCells[faceI];

                    vector NfPatchI = SfPatch[faceI]/magSfPatch[faceI];

                    // Correction vectors for uncoupled boundaries are
                    //  set to zero. We will avoid their usage here.
                    vector DfPatchI = (CfPatch[faceI] - C[own])
                                    * deltaCoeffsPatch[faceI];

                    secAlphaPatch[faceI] = 1.0/(NfPatchI & DfPatchI);
                }
            }
        }
    }

    return tsecAlpha;
}


template<class Type>
Foam::tmp<Foam::surfaceScalarField> Foam::fluxConservative<Type>::weights
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<surfaceScalarField> tw
    (
        new surfaceScalarField
        (
            IOobject
            (
                "fluxConservativeWeightingFactors" + vf.name(),
                mesh.time().timeName(),
                mesh
            ),
            mesh ,
            dimless
        )
    );

    surfaceScalarField& w = tw();
    scalarField& wIn = w.internalField();

    // Mesh addressing
    const unallocLabelList& owner = mesh.owner();

    // Mesh and basic surface interpolation data
    const surfaceScalarField& weights = mesh.weights();
    const scalarField& weightsIn = weights.internalField();

    // Gamma
    const volScalarField& gamma = gamma_;
    const scalarField& gammaIn = gamma.internalField();

    // Linear (!!!) interpolated gamma
    tmp<surfaceScalarField> tgammaf(linearInterpolate(gamma));
    const surfaceScalarField& gammaf = tgammaf();
    const scalarField& gammafIn = gammaf.internalField();

    // Calculate internal w
    forAll (owner, faceI)
    {
        label own = owner[faceI];

        wIn[faceI] = gammaIn[own]/gammafIn[faceI] * weightsIn[faceI];
    }

    // Calculate boundary w
    forAll (mesh.boundary(), patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];

        scalarField& wPatch = w.boundaryField()[patchI];

        // Coupled patches
        if (patch.coupled())
        {
            const unallocLabelList& faceCells = patch.patch().faceCells();

            const scalarField& weightsPatch = weights.boundaryField()[patchI];

            const scalarField& gammafPatch = gammaf.boundaryField()[patchI];

            forAll (patch, faceI)
            {
                const label own = faceCells[faceI];

                wPatch[faceI] =
                    gammaIn[own]/gammafPatch[faceI] * weightsPatch[faceI];
            }
        }
        else
        {
            // Boundary w for uncoupled patches are 1
            wPatch = 1;
        }
    }

    return tw;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> >
Foam::fluxConservative<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    if (!jumpFluxPtr_)
    {
        FatalErrorIn
        (
            "tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > fvmGrad\n"
            "Foam::fluxConservative<Type>::correction\n"
            "(\n"
            "    GeometricField<Type, fvPatchField, volMesh>&"
            ")\n"
        )   << "Jump flux pointer not assigned.\n"
            << "Maybe this scheme was constructed from IStream without\n"
            << " getting a name for the jump flux field?"
            << " If this is true, ::correction(vf) must not be used"
            << " and ::corrected() should actually return 'false'!"
            << abort(FatalError);
    }

    const fvMesh& mesh = this->mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tcorr
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "fluxConservativeCorrection" + vf.name(),
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimensioned<Type>(word(), vf.dimensions(), pTraits<Type>::zero)
        )
    );

    GeometricField<Type, fvsPatchField, surfaceMesh>& corr = tcorr();
    Field<Type>& corrIn = corr.internalField();

    // Mesh addressing
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    // Mesh and basic surface interpolation data
    const surfaceScalarField& weights = mesh.weights();
    const scalarField& weightsIn = weights.internalField();
    const surfaceScalarField& deltaCoeffs = mesh.deltaCoeffs();
    const scalarField& deltaCoeffsIn = deltaCoeffs.internalField();
    const surfaceScalarField& magSf = mesh.magSf();
    const scalarField& magSfIn = magSf.internalField();

    // Non-orthogonal correction secant factors for one-sided gradient.
    const surfaceScalarField secAlpha = secAlphaNonOrth()();
    const scalarField& secAlphaIn = secAlpha.internalField();

    // Gamma
    const volScalarField& gamma = gamma_;
    const scalarField& gammaIn = gamma.internalField();

    // Linear (!!!) interpolated gamma
    tmp<surfaceScalarField> tgammaf(linearInterpolate(gamma));
    const surfaceScalarField& gammaf = tgammaf();
    const scalarField& gammafIn = gammaf.internalField();

    // Jump flux
    const GeometricField<Type, fvsPatchField, surfaceMesh>& jumpFlux =
        *jumpFluxPtr_;
    const Field<Type>& jumpFluxIn = jumpFlux.internalField();

    // Calculate internal correction
    forAll (owner, faceI)
    {
        label own = owner[faceI];
        label nei = neighbour[faceI];

        scalar gammaRelDiff = (gammaIn[nei] - gammaIn[own])/gammafIn[faceI];

        scalar wPwN = weightsIn[faceI] * (1 - weightsIn[faceI]);

        scalar dByMagSf = 1.0/deltaCoeffsIn[faceI]/magSfIn[faceI];

        corrIn[faceI] = secAlphaIn[faceI]
                      * gammaRelDiff * wPwN * dByMagSf
                      * jumpFluxIn[faceI];
    }

    // Calculate boundary correction
    // (Correction will be zero if grad(gamma)*n == 0)
    forAll (mesh.boundary(), patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];

        Field<Type>& corrPatch = corr.boundaryField()[patchI];

        const unallocLabelList& faceCells = patch.patch().faceCells();

        const scalarField& weightsPatch = weights.boundaryField()[patchI];
        const scalarField& deltaCoeffsPatch = deltaCoeffs.boundaryField()[patchI];
        const scalarField& magSfPatch = magSf.boundaryField()[patchI];

        const scalarField& secAlphaPatch = secAlpha.boundaryField()[patchI];

        const scalarField& gammaPatch = gamma.boundaryField()[patchI];
        const scalarField& gammafPatch = gammaf.boundaryField()[patchI];
        const Field<Type>& jumpFluxPatch = jumpFlux.boundaryField()[patchI];

        forAll (patch, faceI)
        {
            const label own = faceCells[faceI];

            scalar gammaRelDiff =
                (gammaPatch[faceI] - gammaIn[own])/gammafPatch[faceI];

            scalar wPwN = weightsPatch[faceI] * (1 - weightsPatch[faceI]);

            scalar dByMagSf = 1.0/deltaCoeffsPatch[faceI]/magSfPatch[faceI];

            corrPatch[faceI] = secAlphaPatch[faceI]
                             * gammaRelDiff * wPwN * dByMagSf
                             * jumpFluxPatch[faceI];
        }
    }

    return tcorr;
}


// ************************************************************************* //
