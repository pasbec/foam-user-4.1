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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
const IOobject
regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::
regionIO
(
    const label& regionI,
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
regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::
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
    size_(regionMesh.size()),
    regionNames_(regionMesh.regionNames()),
    regionMesh_(regionMesh),
    fields_(List<GeometricField<Type, PatchField, GeoMesh>*>(size_,NULL))
{
    forAll (regionNames_, regionI)
    {
        if (debug)
        {
            Info << "regionGeometricField::regionGeometricField(...) : "
                << "Create field for region "
                << regionMesh_.regionName(regionI)
                << endl;
        }

        fields_[regionI] =
            new GeometricField<Type, PatchField, GeoMesh>
            (
                regionIO(regionI, io),
                regionMesh_[regionI],
                dim,
                patchFieldType
            );
    }
}

template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::
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
    size_(regionMesh.size()),
    regionNames_(regionMesh.regionNames()),
    regionMesh_(regionMesh),
    fields_(List<GeometricField<Type, PatchField, GeoMesh>*>(size_,NULL))
{
    forAll (regionNames_, regionI)
    {
        if (debug)
        {
            Info << "regionGeometricField::regionGeometricField(...) : "
                << "Create field for region "
                << regionMesh_.regionName(regionI)
                << endl;
        }

        fields_[regionI] =
            new GeometricField<Type, PatchField, GeoMesh>
            (
                regionIO(regionI, io),
                regionMesh_[regionI],
                dim,
                patchFieldTypes
            );
    }
}

template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::
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
    size_(regionMesh.size()),
    regionNames_(regionMesh.regionNames()),
    regionMesh_(regionMesh),
    fields_(List<GeometricField<Type, PatchField, GeoMesh>*>(size_,NULL))
{
    forAll (regionNames_, regionI)
    {
        if (debug)
        {
            Info << "regionGeometricField::regionGeometricField(...) : "
                << "Create field for region "
                << regionMesh_.regionName(regionI)
                << endl;
        }

        fields_[regionI] =
            new GeometricField<Type, PatchField, GeoMesh>
            (
                regionIO(regionI, io),
                regionMesh_[regionI]
            );
    }
}

template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::
regionGeometricField
(
    const IOobject& io,
    const regionGeometricField
    <
        Type, PatchField, GeoMesh,
        RegionGeoMesh
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
    size_(rgf.mesh().size()),
    regionNames_(rgf.mesh().regionNames()),
    regionMesh_(rgf.mesh()),
    fields_(List<GeometricField<Type, PatchField, GeoMesh>*>(size_,NULL))
{
    forAll (regionNames_, regionI)
    {
        if (debug)
        {
            Info << "regionGeometricField::regionGeometricField(...) : "
                << "Create field for region "
                << regionMesh_.regionName(regionI)
                << endl;
        }

        fields_[regionI] =
            new GeometricField<Type, PatchField, GeoMesh>
            (
                regionIO(regionI, io),
                rgf.field(regionI)
            );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::
~regionGeometricField()
{
    forAll (regionNames_, regionI)
    {
        if (debug)
        {
            Info << "regionGeometricField::~regionGeometricField() : "
                << "Delete field for region "
                << regionMesh_.regionName(regionI)
                << endl;
        }

        delete fields_[regionI];
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
void
regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::operator=
(
    const regionGeometricField
        <
            Type,
            PatchField,
            GeoMesh,
            RegionGeoMesh
        >& rgf
)
{
    if (this == &rgf)
    {
        FatalErrorIn
        (
            "regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::operator="
            "(const regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>&)"
        )   << "attempted assignment to self"
            << abort(FatalError);
    }

    checkField(*this, rgf, "=");

    forAll (this->mesh().regionNames(), regionI)
    {
        checkField(this->field(regionI), rgf.field(regionI), "=");

        // only equate field contents not ID

        this->field(regionI).dimensionedInternalField()
            = rgf.field(regionI).dimensionedInternalField();

        this->field(regionI).boundaryField()
            = rgf.field(regionI).boundaryField();
    }
}

template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
void
regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::operator=
(
    const tmp
    <
        regionGeometricField
        <
            Type, PatchField, GeoMesh,
            RegionGeoMesh
        >
    >& trgf
)
{
    if (this == &(trgf()))
    {
        FatalErrorIn
        (
            "regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::operator="
            "(const tmp<regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh> >&)"
        )   << "attempted assignment to self"
            << abort(FatalError);
    }

    const regionGeometricField
    <
        Type, PatchField, GeoMesh,
        RegionGeoMesh
    >&
    rgf = trgf();

    checkField(*this, rgf, "=");

    forAll (this->mesh().regionNames(), regionI)
    {
        checkField(this->field(regionI), rgf.field(regionI), "=");

        // only equate field contents not ID

        this->field(regionI).dimensions()
            = rgf.field(regionI).dimensions();

        this->field(regionI).internalField().transfer
        (
            const_cast<Field<Type>&>(rgf.field(regionI).internalField())
        );

        this->field(regionI).boundaryField()
            = rgf.field(regionI).boundaryField();
    }

    trgf.clear();
}

// TODO: Parallel?
template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
void
regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::mapBoundaryField
(
    const label& regionI
) const
{
    const GeometricField<Type, PatchField, GeoMesh>& vf0 =
        field(polyMesh::defaultRegion);

    GeometricField<Type, PatchField, GeoMesh>& vf = field(regionI);

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

                const labelList& fmap = mesh().faceMap(regionI);

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
};

// TODO: Parallel?
template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
void
regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::copyInternalBoundaryField
(
    const label& regionI
) const
{
    GeometricField<Type, PatchField, GeoMesh>& vf = field(regionI);
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
};

// TODO: Parallel?
template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
void
regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::interpolateBoundaryField
(
    const label& regionI
) const
{
    const GeometricField<Type, PatchField, GeoMesh>& vf0 =
        field(polyMesh::defaultRegion);

    GeometricField<Type, PatchField, GeoMesh>& vf = field(regionI);

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

                const labelList& fmap = mesh().faceMap(regionI);

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
};

template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
void
regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::extrapolateBoundaryField
(
    const label& regionI
) const
{
    GeometricField<Type, PatchField, GeoMesh>& vf = field(regionI);

    fvc::extrapolate(vf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
