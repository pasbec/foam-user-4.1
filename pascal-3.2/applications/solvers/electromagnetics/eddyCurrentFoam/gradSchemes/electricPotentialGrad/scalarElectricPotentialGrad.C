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

    // Mesh references
    const unallocLabelList& own = mesh.owner();
    const unallocLabelList& nei = mesh.neighbour();

    const vectorField& Sf = mesh.Sf();
    const scalarField& magSf = mesh.magSf();
    const scalarField& w = mesh.weights();
    const scalarField& d = mesh.deltaCoeffs();
    const scalarField& V = mesh.V();

    // Derived mesh data
    const vectorField Sn = Sf/magSf;

    const scalarField& vsfIn = vsf.internalField();

    // Prepare gradient data
    tmp<volVectorField> tepGrad
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
    volVectorField& epGrad = tepGrad();
    vectorField& epGradIn = epGrad.internalField();

    // Get list of conductor faces
    labelList cnFaces
    (
        mesh.lookupObject<faceSet>
        (
            nameInterface_
        ).toc()
    );

    // Remove wrongly interpolated gradient
    // contributions from cells at the interface
    // between non-conducting and conducting region
    if (this->baseGradScheme_ == this->GAUSS)
    {
        tmp<surfaceScalarField> tssf
        (
            tinterpScheme_().interpolate(vsf)
        );
        surfaceScalarField& ssf = tssf();
        scalarField& ssfIn = ssf.internalField();

        // Copy Gauss gradient values
        tmp<volVectorField> tgGrad
        (
            gaussGrad<scalar>::gradf
            (
                tssf,
                name
            )
        );
        tepGrad() == tgGrad();

        register label faceI, ownFaceI, neiFaceI;

        // Remove wrong face contributions
        // from interal faces
        forAll(cnFaces, facei)
        {
            faceI = cnFaces[facei];

            ownFaceI = own[faceI];
            neiFaceI = nei[faceI];

            {
                vector Sfssf = Sf[faceI] * ssfIn[faceI];

                epGradIn[ownFaceI] -= Sfssf / V[ownFaceI];
                epGradIn[neiFaceI] += Sfssf / V[neiFaceI];
            }
        }

        // TODO: Remove wrong contributions
        // from boundary faces

        tgGrad.clear();
        tssf.clear();
    }
//     else if (this->baseGradScheme_ == this->LEASTSQUARES)
//     {
//         tmp<leastSquaresGrad<scalar> > tlsScheme
//         (
//             new leastSquaresGrad<scalar>(mesh)
//         );
//         leastSquaresGrad<scalar>& lsScheme = tlsScheme();
//
//         // Copy leastSquares gradient values
//         tmp<volVectorField> tlsGrad
//         (
//             lsScheme.calcGrad
//             (
//                 vsf,
//                 name
//             )
//         );
//         tepGrad() == tlsGrad();
//
//         const leastSquaresVectors& lsv = leastSquaresVectors::New(mesh);
//
//         const surfaceVectorField& ownLs = lsv.pVectors();
//         const surfaceVectorField& neiLs = lsv.nVectors();
//
//         const vectorField& ownLsIn = ownLs.internalField();
//         const vectorField& neiLsIn = neiLs.internalField();
//
//         register label faceI, ownFaceI, neiFaceI;
//
//         // Remove wrong contributions
//         // from interal faces
//         forAll(cnFaces, facei)
//         {
//             faceI = cnFaces[facei];
//
//             ownFaceI = own[faceI];
//             neiFaceI = nei[faceI];
//
//             {
//                 scalar deltaVsf = vsfIn[neiFaceI] - vsfIn[ownFaceI];
//
//                 epGradIn[ownFaceI] -= ownLsIn[faceI] * deltaVsf;
//                 epGradIn[neiFaceI] += neiLsIn[faceI] * deltaVsf;
//             }
//         }
//
//         // TODO: Remove wrong contributions
//         // from boundary faces
//
//         tlsGrad.clear();
//         tlsScheme.clear();
//     }

    // Add correct face value contributions to cells
    // for the electric scalar potential from its
    // corresponding gradient at the interface between
    // non-conducting and conducting region depending on
    // the magnetic vector potential
    {
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

        register label faceI, ownFaceI, neiFaceI;

        forAll(cnFaces, facei)
        {
            faceI = cnFaces[facei];

            ownFaceI = own[faceI];
            neiFaceI = nei[faceI];

            {
                // Face centre to cell centre distances
                scalar fN = w[faceI] / d[faceI];
                scalar Pf = 1.0 / d[faceI] - fN;

                // Gradient weight coeffs
                scalar wP = omega.value() * w[faceI];
                scalar wN = omega.value() - wP;

                // TODO: Interpolation of A?
                // Currently, a linear interpolation is hard-coded!

                // Weighted face gradients
                scalar SfAwOwn = wP * As * (A[ownFaceI] & Sn[faceI]);
                scalar SfAwNei = wN * As * (A[neiFaceI] & Sn[faceI]);

                // Extrapolated face values
                scalar ssfOwn = vsf[ownFaceI] + Pf * (SfAwOwn + SfAwNei);
                scalar ssfNei = vsf[neiFaceI] - fN * (SfAwOwn + SfAwNei);

                epGrad[ownFaceI] += Sf[faceI] * ssfOwn / V[ownFaceI];
                epGrad[neiFaceI] -= Sf[faceI] * ssfNei / V[neiFaceI];
            }
        }
    }

    epGrad.rename("grad(" + vsf.name() + ')');
    this->correctBoundaryConditions(vsf, epGrad);

    return tepGrad;
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
