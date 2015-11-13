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

#include "regionGeometricField.H"

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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

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
    forAll(regionNames_, regionI)
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
    forAll(regionNames_, regionI)
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
    forAll(regionNames_, regionI)
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
    forAll(regionNames_, regionI)
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

    forAll(this->mesh().regionNames(), regionI)
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

    forAll(this->mesh().regionNames(), regionI)
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
