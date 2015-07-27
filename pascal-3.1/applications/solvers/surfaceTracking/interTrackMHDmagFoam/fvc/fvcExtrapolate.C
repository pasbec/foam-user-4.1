/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

#include "beGaussGrad.H"

#include "emptyFvPatchFields.H"
#include "symmetryFvPatchFields.H"
#include "wedgeFvPatchFields.H"
#include "calculatedFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"

// TODO [High] One could also realize this in boundary condition?

// TODO [High] Currently only for calculated patches and fails if any other patch type is use
//             except for "special" patches - which are ignored anyway. There might be a better
//             approach?

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void
extrapolate
(
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const scalar& alpha
)
{
    if (alpha < 0.0 || alpha > 1.0)
    {
        FatalErrorIn("fvc::extrapolate(...)")
            << "Value for alpha out of range. The weighting factor "
            << " alpha may only have values between 0.0 and 1.0."
            << abort(FatalError);
    }

    const fvMesh& mesh = vf.mesh();

    typedef typename outerProduct<vector, Type>::type GradType;

    // Create copy of vf with zero-gradient patches
    tmp< GeometricField<Type, fvPatchField, volMesh> > tvfZg
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                vf.name(),
                mesh.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            vf,
            zeroGradientFvPatchField<Type>::typeName
        )
    );

    GeometricField<Type, fvPatchField, volMesh>& vfZg = tvfZg();

    // Create boundary extrapolated gradient scheme
    tmp<fv::beGaussGrad<Type> > tvfGradScheme
    (
        new fv::beGaussGrad<Type>(mesh)
     );

    fv::beGaussGrad<Type>& vfGradScheme = tvfGradScheme();

    // Gradient field of vf with extrapolated boundary values
    tmp<GeometricField<GradType, fvPatchField, volMesh> > tvfGrad
    (
        new GeometricField<GradType, fvPatchField, volMesh>
        (
            vfGradScheme.grad(vfZg)
        )
    );

    GeometricField<GradType, fvPatchField, volMesh>& vfGrad = tvfGrad();

    tvfZg.clear();
    tvfGradScheme.clear();

    // Extrapolate boundary values
    forAll(vf.boundaryField(), patchI)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchI];

        if
        (
            vf.boundaryField()[patchI].type()
            == calculatedFvPatchField<Type>::typeName
        )
        {
            calculatedFvPatchField<Type>& vfp =
                refCast<calculatedFvPatchField<Type> >
                (
                    vf.boundaryField()[patchI]
                );

            tmp<Field<vector> > deltap =
                mesh.boundary()[patchI].delta();

            tmp<Field<Type> > vfpi =
                vf.boundaryField()[patchI].patchInternalField();

            tmp<Field<GradType> > vfpGrad = alpha
              * vfGrad.boundaryField()[patchI];

            tmp<Field<GradType> > vfpiGrad = (1.0-alpha)
              * vfGrad.boundaryField()[patchI].patchInternalField();

            // Extrapolate boundary values
            vfp == vfpi + ((vfpGrad+vfpiGrad) & deltap);

            vfpiGrad.clear();
            vfpGrad.clear();
            vfpi.clear();
            deltap.clear();
        }
        else if
        (
            (pp.type() != emptyFvPatch::typeName)
         && (pp.type() != symmetryFvPatch::typeName)
         && (pp.type() != wedgeFvPatch::typeName)
         && !pp.coupled()
        )
        {
            FatalErrorIn("fvc::extrapolate(...)")
                << "Boundary patch "
                << pp.name()
                << " of field " << vf.name()
                << " is not of type \"calculated\"."
                << abort(FatalError);
        }
    }

    tvfGrad.clear();
}

template<class Type>
void
extrapolate
(
    tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
    const scalar& alpha
)
{
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf();

    extrapolate(vf);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
