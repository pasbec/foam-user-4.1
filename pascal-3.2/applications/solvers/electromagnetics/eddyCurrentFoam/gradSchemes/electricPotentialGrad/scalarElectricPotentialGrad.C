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

    // Derived mesh data
    const vectorField Sn = Sf/magSf;

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

    //
    if (this->baseScheme() == this->GAUSS)
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

        // TODO: Check if vsf.name() ~ nameElectricPotential_

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

        const volScalarField& sigma =
            mesh.lookupObject<volScalarField>
            (
                nameConductivity_
            );

        const volVectorField& A =
            mesh.lookupObject<volVectorField>
            (
                nameMagneticPotential_ + complexPartA
            );

        scalar As = complexSignA;

        // Internal contributions
        forAll(own, faceI)
        // TODO: Optimize speed
        {
            label ownFaceI = own[faceI];
            label neiFaceI = nei[faceI];

            {
                // Interpolation weights
                scalar wP = w[faceI];
                scalar wN = 1.0 - wP;

                // Conductivity
                scalar sigmaOwn = sigma[ownFaceI];
                scalar SigmaNei = sigma[neiFaceI];

                // Conductivity-weighted inverse face centre
                // to cell centre distances
                scalar dByWsigmaOwn = d[faceI] / wN * sigmaOwn;
                scalar dByWsigmaNei = d[faceI] / wP * SigmaNei;

                // Compex signed sigma difference
                scalar sigmaDiff = As * (sigmaOwn - SigmaNei);

                // Face interpolated complex time derivative of A
                scalar wPownAn = wP * (A[ownFaceI] & Sn[faceI]);
                scalar wNneiAn = wN * (A[neiFaceI] & Sn[faceI]);
                scalar ddtAfn = omega.value() * (wPownAn + wNneiAn);

                // Current conserving face value of V
                scalar Vf = sigmaDiff * ddtAfn;
                Vf += dByWsigmaOwn * vsf[ownFaceI];
                Vf += dByWsigmaNei * vsf[neiFaceI];
                Vf /= dByWsigmaOwn + dByWsigmaNei;

                // Gradient contributions
                epGrad[ownFaceI] += Sf[faceI] * Vf;
                epGrad[neiFaceI] -= Sf[faceI] * Vf;
            }
        }

        // Boundary contributions
        // TEST
        // TODO: Is this enough? What if there are symmetry,
        //       empty or other special bc?
        {
            tmp<surfaceScalarField> tssf
            (
                tinterpScheme_().interpolate(vsf)
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
                    epGrad[pFaceCells[facei]] += pSf[facei]*pssf[facei];
                }
            }

            tssf.clear();
        }

        epGradIn /= mesh.V();
    }
    else if (this->baseScheme() == this->LEASTSQUARES)
    {
        FatalErrorIn
        (
            "electricPotentialGrad<scalar>::calcGrad\n"
            "(\n"
            "    GeometricField<Type, fvPatchField, volMesh>&\n"
            ")\n"
        )   << "Gradient calculation for LEASTSQUARES "
            << "base scheme not yet implemented."
            << abort(FatalError);
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
    const GeometricField<scalar, fvPatchField, volMesh>& vsf
) const
{
    FatalErrorIn
    (
        "electricPotentialGrad<scalar>::fvmGrad\n"
        "(\n"
        "    GeometricField<Type, fvPatchField, volMesh>&\n"
        ")\n"
    )   << "Implicit gradient operator not yet implemented."
        << abort(FatalError);

    tmp<BlockLduSystem<vector, vector> > tbs
    (
        new BlockLduSystem<vector, vector>(vsf.mesh())
    );
    return tbs;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
