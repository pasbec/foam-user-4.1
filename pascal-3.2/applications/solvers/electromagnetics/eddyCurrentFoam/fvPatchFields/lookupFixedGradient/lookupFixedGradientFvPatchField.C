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

#include "lookupFixedGradientFvPatchField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
lookupFixedGradientFvPatchField<Type>::lookupFixedGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    gradient_(p.size(), pTraits<Type>::zero),
    lookupName_("none"),
    lookupType_("surface")
{}


template<class Type>
lookupFixedGradientFvPatchField<Type>::lookupFixedGradientFvPatchField
(
    const lookupFixedGradientFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
    gradient_(ptf.gradient_, mapper),
    lookupName_(ptf.lookupName_),
    lookupType_(ptf.lookupType_)
{}


template<class Type>
lookupFixedGradientFvPatchField<Type>::lookupFixedGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict),
    gradient_("gradient", dict, p.size()),
    lookupName_(dict.lookupOrDefault<word>("lookupName", word())),
    lookupType_(dict.lookupOrDefault<word>("lookupType", "surface"))
{
    evaluate();
}


template<class Type>
lookupFixedGradientFvPatchField<Type>::lookupFixedGradientFvPatchField
(
    const lookupFixedGradientFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),
    gradient_(ptf.gradient_),
    lookupName_(ptf.lookupName_),
    lookupType_(ptf.lookupType_)
{}


template<class Type>
lookupFixedGradientFvPatchField<Type>::lookupFixedGradientFvPatchField
(
    const lookupFixedGradientFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    gradient_(ptf.gradient_),
    lookupName_(ptf.lookupName_),
    lookupType_(ptf.lookupType_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void lookupFixedGradientFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<Type>::autoMap(m);
    gradient_.autoMap(m);
}


template<class Type>
void lookupFixedGradientFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    const lookupFixedGradientFvPatchField<Type>& fgptf =
        refCast<const lookupFixedGradientFvPatchField<Type> >(ptf);

    gradient_.rmap(fgptf.gradient_, addr);
}


template<class Type>
void lookupFixedGradientFvPatchField<Type>::updateCoeffs()
{
}


template<class Type>
void lookupFixedGradientFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type>::operator=
    (
        this->patchInternalField() + gradient_/this->patch().deltaCoeffs()
    );

    fvPatchField<Type>::evaluate();
}


template<class Type>
tmp<Field<Type> > lookupFixedGradientFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::one)
    );
}


template<class Type>
tmp<Field<Type> > lookupFixedGradientFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return gradient()/this->patch().deltaCoeffs();
}


template<class Type>
tmp<Field<Type> > lookupFixedGradientFvPatchField<Type>::
gradientInternalCoeffs() const
{
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}


template<class Type>
tmp<Field<Type> > lookupFixedGradientFvPatchField<Type>::
gradientBoundaryCoeffs() const
{
    return gradient();
}


template<class Type>
void lookupFixedGradientFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    gradient_.writeEntry("gradient", os);
    os.writeKeyword("lookupName") << lookupName_ << token::END_STATEMENT << nl;
    os.writeKeyword("lookupType") << lookupType_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
