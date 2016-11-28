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

InNamespace
    Foam::fvc

Description
    Surface integrate two surfaceFields creating a volField.
    Surface sum two surfaceFields creating a volField.

SourceFiles
    fvcJumpSurfaceIntegrate.C

\*---------------------------------------------------------------------------*/


#ifndef fvcJumpSurfaceIntegrate_H
#define fvcJumpSurfaceIntegrate_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
    template<class Type>
    void jumpSurfaceIntegrate
    (
        Field<Type>&,
        const GeometricField<Type, fvsPatchField, surfaceMesh>&,
        const GeometricField<Type, fvsPatchField, surfaceMesh>&
    );

    template<class Type>
    tmp<GeometricField<Type, fvPatchField, volMesh> >
    jumpSurfaceIntegrate
    (
        const GeometricField<Type, fvsPatchField, surfaceMesh>&,
        const GeometricField<Type, fvsPatchField, surfaceMesh>&
    );

    template<class Type>
    tmp<GeometricField<Type, fvPatchField, volMesh> >
    jumpSurfaceIntegrate
    (
        const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >&,
        const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >&
    );

    template<class Type>
    tmp<GeometricField<Type, fvPatchField, volMesh> > jumpSurfaceSum
    (
        const GeometricField<Type, fvsPatchField, surfaceMesh>&,
        const GeometricField<Type, fvsPatchField, surfaceMesh>&
    );

    template<class Type>
    tmp<GeometricField<Type, fvPatchField, volMesh> > jumpSurfaceSum
    (
        const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >&,
        const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fvcJumpSurfaceIntegrate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
