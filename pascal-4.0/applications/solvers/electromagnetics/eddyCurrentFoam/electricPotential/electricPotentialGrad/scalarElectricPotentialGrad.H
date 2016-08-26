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

Typedef
    Foam::fv::scalarElectricPotentialGrad

Description
    Specialisation of electricPotentialGrad for scalars.

SourceFiles
    scalarElectricPotentialGrad.C

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#ifndef scalarElectricPotentialGrad_H
#define scalarElectricPotentialGrad_H

#include "electricPotentialGrad.H"

// TODO: Realize this with simple interpolation scheme instead of a
//       dedicated gradient scheme?
//       Is it possible to implement an implicit block-matrix version in that
//       case?

// TODO: Implicit gradient

// TODO: Check delta coeffs and consistency with electricPotentialLaplacian
//

// TODO: Check boundary contributions for boundaries with transformation!!
//       Currently boundary contributions are simply realized to work with
//       coupled processor patches

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
) const;


template<>
tmp<BlockLduSystem<vector, vector> > electricPotentialGrad<scalar>::fvmGrad
(
    const GeometricField<scalar, fvPatchField, volMesh>&
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
