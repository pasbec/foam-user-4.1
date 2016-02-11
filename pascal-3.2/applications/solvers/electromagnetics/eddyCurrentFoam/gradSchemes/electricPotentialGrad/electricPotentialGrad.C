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

\*---------------------------------------------------------------------------*/

#include "electricPotentialGrad.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void electricPotentialGrad<Type>::readMagneticProperties()
{
    dictionary namesDict(magneticProperties_.subDict("names"));

    nameInterface_ = word(namesDict.lookup("interface"));
    nameConductivity_ = word(namesDict.lookup("conductivity"));
    nameElectricPotential_ = word(namesDict.lookup("electricPotential"));
    nameMagneticPotential_ = word(namesDict.lookup("magneticPotential"));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >
>
electricPotentialGrad<Type>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    FatalErrorIn
    (
        "tmp<GeometricField> calcGrad\n"
        "(\n"
        "    GeometricField<Type, fvPatchField, volMesh>&,\n"
        "    const word&,\n"
        ")\n"
    )   << "Gradient calculation defined only for scalar."
        << abort(FatalError);

    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = vsf.mesh();

    tmp<GeometricField<GradType, fvPatchField, volMesh> > tepGrad
    (
        new GeometricField<GradType, fvPatchField, volMesh>
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
            dimensioned<GradType>
            (
                "zero",
                vsf.dimensions()/dimLength,
                pTraits<GradType>::zero
            ),
            zeroGradientFvPatchField<GradType>::typeName
        )
    );

    return tepGrad;
}


template<class Type>
tmp
<
    BlockLduSystem<vector, typename outerProduct<vector, Type>::type>
>
electricPotentialGrad<Type>::fvmGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    FatalErrorIn
    (
        "tmp<BlockLduSystem> fvmGrad\n"
        "(\n"
        "    GeometricField<Type, fvPatchField, volMesh>&\n"
        ")\n"
    )   << "Implicit gradient operator defined only for scalar."
        << abort(FatalError);

    typedef typename outerProduct<vector, Type>::type GradType;

    tmp<BlockLduSystem<vector, GradType> > tbs
    (
        new BlockLduSystem<vector, GradType>(vf.mesh())
    );

    return tbs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
