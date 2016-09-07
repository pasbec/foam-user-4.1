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

#include "regionGeometricField.H"
#include "fvcExtrapolate.H"
#include "calculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// check mesh for two fields
#define checkField(gf1, gf2, op)                                    \
if ((gf1).mesh() != (gf2).mesh())                                   \
{                                                                   \
    FatalErrorIn("checkField(gf1, gf2, op)")                        \
        << "different mesh for fields "                             \
        << (gf1).name() << " and " << (gf2).name()                  \
        << " during operatrion " <<  op                             \
        << abort(FatalError);                                       \
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
const Foam::IOobject
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::
regionIO
(
    label regionI,
    const IOobject& io
) const
{
    const objectRegistry& db = regionMesh_[regionI].thisDb();

    return IOobject
    (
        io.name(),
        io.time().timeName(),
        db,
        io.readOpt(),
        io.writeOpt()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::
regionGeometricField
(
    const IOobject& io,
    const RegionMesh& regionMesh,
    const dimensioned<Type>& dim,
    const word& patchFieldType
)
:
    regIOobject
    (
        IOobject
        (
            io.name(),
            io.time().timeName(),
            regionMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        )
    ),
    regionMesh_(regionMesh),
    fieldPtrs_(regionMesh.size())
{
    forAll (*this, regionI)
    {
        if (debug)
        {
            Info<< "Foam::regionGeometricField::regionGeometricField(...) : "
                << "Create field for region "
                << regionMesh_.regions()[regionI]
                << endl;
        }

        fieldPtrs_.set
        (
            regionI,
            new GeometricField<Type, PatchField, GeoMesh>
            (
                regionIO(regionI, io),
                regionMesh_[regionI],
                dim,
                patchFieldType
            )
        );
    }
}


template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::
regionGeometricField
(
    const IOobject& io,
    const RegionMesh& regionMesh,
    const dimensioned<Type>& dim,
    const wordList& patchFieldTypes
)
:
    regIOobject
    (
        IOobject
        (
            io.name(),
            io.time().timeName(),
            regionMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        )
    ),
    regionMesh_(regionMesh),
    fieldPtrs_(regionMesh.size())
{
    forAll (*this, regionI)
    {
        if (debug)
        {
            Info<< "Foam::regionGeometricField::regionGeometricField(...) : "
                << "Create field for region "
                << regionMesh_.regions()[regionI]
                << endl;
        }

        fieldPtrs_.set
        (
            regionI,
            new GeometricField<Type, PatchField, GeoMesh>
            (
                regionIO(regionI, io),
                regionMesh_[regionI],
                dim,
                patchFieldTypes
            )
        );
    }
}


template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::
regionGeometricField
(
    const IOobject& io,
    const RegionMesh& regionMesh
)
:
    regIOobject
    (
        IOobject
        (
            io.name(),
            io.time().timeName(),
            regionMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        )
    ),
    regionMesh_(regionMesh),
    fieldPtrs_(regionMesh.size())
{
    forAll (*this, regionI)
    {
        if (debug)
        {
            Info<< "Foam::regionGeometricField::regionGeometricField(...) : "
                << "Create field for region "
                << regionMesh_.regions()[regionI]
                << endl;
        }

        fieldPtrs_.set
        (
            regionI,
            new GeometricField<Type, PatchField, GeoMesh>
            (
                regionIO(regionI, io),
                regionMesh_[regionI]
            )
        );
    }
}


template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::
regionGeometricField
(
    const IOobject& io,
    const regionGeometricField
    <
        Type, PatchField, GeoMesh, RegionGeoMesh
    >& rgf
)
:
    regIOobject
    (
        IOobject
        (
            io.name(),
            io.time().timeName(),
            rgf.mesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        )
    ),
    regionMesh_(rgf.mesh()),
    fieldPtrs_(rgf.mesh().size())
{
    forAll (*this, regionI)
    {
        if (debug)
        {
            Info<< "Foam::regionGeometricField<>::regionGeometricField(...) : "
                << "Create field for region "
                << regionMesh_.regions()[regionI]
                << endl;
        }

        fieldPtrs_.set
        (
            regionI,
            new GeometricField<Type, PatchField, GeoMesh>
            (
                regionIO(regionI, io),
                rgf.field(regionI)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::
~regionGeometricField()
{
    fieldPtrs_.clear();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
void
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::operator=
(
    const regionGeometricField
        <
            Type, PatchField, GeoMesh, RegionGeoMesh
        >& rgf
)
{
    if (this == &rgf)
    {
        FatalErrorIn
        (
            "Foam::regionGeometricField<>::operator=(...)"
        )   << "Attempted assignment to self"
            << abort(FatalError);
    }

    checkField(*this, rgf, "=");

    forAll (*this, regionI)
    {
        checkField(operator[](regionI), rgf.field(regionI), "=");

        // only equate field contents not ID

        operator[](regionI).dimensionedInternalField()
            = rgf.field(regionI).dimensionedInternalField();

        operator[](regionI).boundaryField()
            = rgf.field(regionI).boundaryField();
    }
}


template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
void
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::operator=
(
    const tmp
    <
        regionGeometricField
        <
            Type, PatchField, GeoMesh, RegionGeoMesh
        >
    >& trgf
)
{
    if (this == &(trgf()))
    {
        FatalErrorIn
        (
            "Foam::regionGeometricField<>::operator=(...)"
        )   << "Attempted assignment to self"
            << abort(FatalError);
    }

    const regionGeometricField
    <
        Type, PatchField, GeoMesh,
        RegionGeoMesh
    >&
    rgf = trgf();

    checkField(*this, rgf, "=");

    forAll (*this, regionI)
    {
        checkField(operator[](regionI), rgf.field(regionI), "=");

        // only equate field contents not ID

        operator[](regionI).dimensions()
            = rgf.field(regionI).dimensions();

        operator[](regionI).internalField().transfer
        (
            const_cast<Field<Type>&>(rgf.field(regionI).internalField())
        );

        operator[](regionI).boundaryField()
            = rgf.field(regionI).boundaryField();
    }

    trgf.clear();
}


// * * * * * * * * * * * * * * * * Mapping * * * * * * * * * * * * * * * * * //

// TODO: Parallel?
template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
void
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::
mapBoundaryField
(
    label regionI
) const
{
    label regionI0 = regions()[polyMesh::defaultRegion];

    const GeometricField<Type, PatchField, GeoMesh>& vf0 = operator[](regionI0);

    GeometricField<Type, PatchField, GeoMesh>& vf = operator[](regionI);

    const polyBoundaryMesh& pbm0 = vf0.mesh().boundaryMesh();
    const polyBoundaryMesh& pbm = vf.mesh().boundaryMesh();

    forAll (pbm, patchI)
    {
        if
        (
            vf.boundaryField()[patchI].type()
         == calculatedFvPatchField<Type>::typeName
        )
        {
            const polyPatch& patch = pbm[patchI];

            label patchI0 = pbm0.findPatchID(patch.name());

            // Patch is present in regionI AND also in
            // default region. Mapping is done based on
            // direct addressing.
            if (patchI0 > -1)
            {
                const polyPatch& patch0 = pbm0[patchI0];
                const Field<Type>& patchField0 = vf0.boundaryField()[patchI0];
                label patchStart0 = patch0.start();
                label patchSize0 = patchField0.size();

                Field<Type>& patchField = vf.boundaryField()[patchI];
                label patchStart = patch.start();

                const labelList& fmap =
                    mesh().typeMap(addressingTypes::FACE, regionI, regionI0);

                forAll (patchField, facei)
                {
                    label faceI = patchStart + facei;
                    label faceI0 = fmap[faceI];

                    label facei0 = faceI0 - patchStart0;

                    if (facei0 > -1 && facei0 < patchSize0)
                    {
                        patchField[facei] = patchField0[facei0];
                    }
                    else
                    {
                        // TODO: Mapped from internal face. Do what?
                        //       one side of an intersecting patch?
                        // NOTE: This should already be avoided during
                        //       construction of regionPolyMesh
                        //       (or the first time faceMap is used) by
                        //       checking all boundary faces?
                    }
                }
            }
        }
    }

    vf.correctBoundaryConditions();
};


template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
void
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::
copyInternalBoundaryField
(
    label regionI
) const
{
    GeometricField<Type, PatchField, GeoMesh>& vf = operator[](regionI);
    const polyBoundaryMesh& pbm = vf.mesh().boundaryMesh();

    forAll (pbm, patchI)
    {
        if
        (
            vf.boundaryField()[patchI].type()
         == calculatedFvPatchField<Type>::typeName
        )
        {
            vf.boundaryField()[patchI] ==
                vf.boundaryField()[patchI].patchInternalField();
        }
    }

    vf.correctBoundaryConditions();
};


// TODO: Parallel?
template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
void
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::
interpolateBoundaryField
(
    label regionI
) const
{
    label regionI0 = regions()[polyMesh::defaultRegion];

    const GeometricField<Type, PatchField, GeoMesh>& vf0 = operator[](regionI0);

    GeometricField<Type, PatchField, GeoMesh>& vf = operator[](regionI);

    const polyBoundaryMesh& pbm0 = vf0.mesh().boundaryMesh();
    const polyBoundaryMesh& pbm = vf.mesh().boundaryMesh();

    forAll (pbm, patchI)
    {
        if
        (
            vf.boundaryField()[patchI].type()
         == calculatedFvPatchField<Type>::typeName
        )
        {
            const polyPatch& patch = pbm[patchI];

            label patchI0 = pbm0.findPatchID(patch.name());

            // Patch is only present in regionI but NOT in
            // default region. Thus, all faces of this patch
            // correspond to interal faces of default region.
            // Linear interpolation is applied to get values
            // for the corresponding patch field.
            if (patchI0 == -1)
            {
                const Field<Type>& vf0In = vf0.internalField();
                const scalarField& w0 = vf0.mesh().weights().internalField();
                const labelList& own0 = vf0.mesh().faceOwner();
                const labelList& ngb0 = vf0.mesh().faceNeighbour();

                Field<Type>& patchField = vf.boundaryField()[patchI];
                label patchStart = patch.start();

                const labelList& fmap =
                    mesh().typeMap(addressingTypes::FACE, regionI, regionI0);

                forAll (patchField, facei)
                {
                    label faceI = patchStart + facei;
                    label faceI0 = fmap[faceI];

                    scalar w0I = w0[faceI0];

                    patchField[facei] =
                        w0I * vf0In[own0[faceI0]]
                    + (1.0 - w0I) * vf0In[ngb0[faceI0]];
                }
            }
            else
            {
                // TODO: What if someone changes the name of
                //       one side of an intersecting patch?
                // NOTE: This should already be avoided during
                //       construction of regionPolyMesh
                //       (or the first time faceMap is used) by
                //       checking all boundary faces?
            }
        }
    }

    vf.correctBoundaryConditions();
};


template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
void
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::
extrapolateBoundaryField
(
    label regionI
) const
{
    GeometricField<Type, PatchField, GeoMesh>& vf = operator[](regionI);

    fvc::extrapolate(vf);
};


// ************************************************************************* //
