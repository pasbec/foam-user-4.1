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

// TODO [High]: Seperate between point and cell mapping
//              depending on vol- or point- fields.

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
void
regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::map() const
{
    forAll(mesh().regionNames(), regionI)
    {
        word regionNameI = mesh().regionName(regionI);

        if (regionNameI != polyMesh::defaultRegion)
        {
            map(regionI);
        }
    }
}

template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
void
regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::map
(
    const label& regionI
) const
{
    const GeometricField<Type, PatchField, GeoMesh>& vf0 =
        field(polyMesh::defaultRegion);

    GeometricField<Type, PatchField, GeoMesh>& vf = field(regionI);

    const labelIOList& map = mesh().cellMap(regionI);

    forAll(field(regionI), cellI)
    {
        vf[cellI] = vf0[map[cellI]];
    }
}

template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
void
regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::map
(
    const word& regionName
) const
{
    map(mesh().regionIndex(regionName));
}

template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
void
regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::rmap() const
{
    forAll(mesh().regionNames(), regionI)
    {
        word regionNameI = mesh().regionName(regionI);

        if (regionNameI != polyMesh::defaultRegion)

        {
            rmap(regionI);
        }
    }
}

template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
void
regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::rmap
(
    const label& regionI
) const
{
    forAll(field(regionI), cellI)
    {
        GeometricField<Type, PatchField, GeoMesh>& vf0 =
            field(polyMesh::defaultRegion);

        const GeometricField<Type, PatchField, GeoMesh>& vf = field(regionI);

        const labelIOList& map = mesh().cellMap(regionI);

        vf0[map[cellI]] = vf[cellI];
    }
}

template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
void
regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::rmap
(
    const word& regionName
) const
{
    rmap(mesh().regionIndex(regionName));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
