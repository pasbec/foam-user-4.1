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

Class
    regionGeometricField

Description
    ...

SourceFiles
    regionGeometricField.C

\*---------------------------------------------------------------------------*/

#ifndef regionGeometricField_H
#define regionGeometricField_H

#include "regIOobject.H"
#include "GeometricFields.H"
#include "dimensionedType.H"

#include "regionGeoMesh.H"
#include "regionPolyMesh.H"

// TODO [Low]: This is a template class. Capitalize first letter!

// TODO [High]: Test parallel mapping

// TODO [High]: Implement template specializations for the currently
//              hard coded membemer functions for boundary field mapping of volFields
//              in regionVolFieldMapping.C

// TODO [Idea]: Make this class demand-driven and allow all fields to be region-fields
// TODO [Idea]: Allow default region to be absent. WE NEED A LOT OF CHECKS!!!

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class regionGeometricField Declaration
\*---------------------------------------------------------------------------*/

template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
class regionGeometricField
:
    public refCount,
    public regIOobject
{
public:

    // Public typedefs

        typedef typename RegionGeoMesh::Mesh RegionMesh;


private:

    // Private data

        //- Reference to regionMesh
        const RegionMesh& regionMesh_;

        //- Field pointers
        mutable PtrList<GeometricField<Type, PatchField, GeoMesh> > fieldPtrs_;

    // Private Member Functions

        // Helper functions for constructor

            //- Return IOobject for region with index regionI
            const IOobject
            regionIO
            (
                label regionI,
                const IOobject& io
            ) const;


public:

    //- Runtime type information
    TypeName("regionGeometricField");

    // Public typedefs

        typedef typename GeoMesh::Mesh Mesh;
        typedef typename GeoMesh::BoundaryMesh BoundaryMesh;

        typedef DimensionedField<Type, GeoMesh> DimensionedInternalField;
        typedef Field<Type> InternalField;
        typedef PatchField<Type> PatchFieldType;

    // Constructors

        //- Constructor given IOobject, mesh, dimensioned<Type>
        //  and patch type.
        regionGeometricField
        (
            const IOobject& io,
            const RegionMesh& regionMesh,
            const dimensioned<Type>& dim,
            const word& patchFieldType
        );

        //- Constructor given IOobject, mesh, dimensioned<Type>
        //  and patch types.
        regionGeometricField
        (
            const IOobject& io,
            const RegionMesh& regionMesh,
            const dimensioned<Type>& dim,
            const wordList& patchFieldTypes
        );

        //- Construct and read given IOobject
        regionGeometricField
        (
            const IOobject& io,
            const RegionMesh& regionMesh
        );

        //- Construct as copy resetting IO parameters
        regionGeometricField
        (
            const IOobject& io,
            const regionGeometricField
            <
                Type, PatchField, GeoMesh,
                RegionGeoMesh
            >& rgf
        );

    // Destructor

        virtual ~regionGeometricField();

    // Member Functions

        // Access

            //- Return const reference to global region count
            inline label size() const
            {
                return regionMesh_.size();
            };

            //- Return const reference to regions
            inline const regionList& regions() const
            {
                return regionMesh_.regions();
            };

            //- Return const reference to region mesh
            inline const RegionMesh& mesh() const
            {
                return regionMesh_;
            };

            //- Return dimensions
            inline const dimensionSet& dimensions() const
            {
                return operator[]
                (
                    regions()[polyMesh::defaultRegion]
                ).dimensions();
            };

            //- Return reference to field
            inline GeometricField<Type, PatchField, GeoMesh>& field
            (
                label regionI
            ) const
            {
                return operator[](regionI);
            };

            //- Return reference to field
            inline GeometricField<Type, PatchField, GeoMesh>& field
            (
                const word& regionName = polyMesh::defaultRegion
            ) const
            {
                return operator[](regions()[regionName]);
            };

            //- WriteData member function required by regIOobject
            inline bool writeData(Ostream&) const
            {
                return false;
            };

        // Internal field data processing

            //- Internal field forward mapping
            void mapInteralField(label regionI) const
            {
                notImplemented(type() + "::mapInteralField()");
            };

            //- Internal field reverse mapping
            void rmapInteralField(label regionI) const
            {
                notImplemented(type() + "::rmapInteralField()");
            };

        // Boundary field data processing

//- TODO: Boundary field mapping where possible
//- TODO: Parallel?
            void mapBoundaryField(label regionI) const
            ;
//             {
//                 notImplemented(type() + "::mapBoundaryField()");
//             };

//- TODO: Boundary field copy from internal field
//- TODO: Parallel?
            void copyInternalBoundaryField(label regionI) const
            ;
//             {
//                 notImplemented(type() + "::copyInternalBoundaryField()");
//             };

//- TODO: Boundary field interpolation
//- TODO: Parallel?
            void interpolateBoundaryField(label regionI) const
            ;
//             {
//                 notImplemented(type() + "::interpolateBoundaryField()");
//             };

//- TODO: Boundary field extrapolation
//- TODO: Parallel?
            void extrapolateBoundaryField(label regionI) const
            ;
//             {
//                 notImplemented(type() + "::extrapolateBoundaryField()");
//             };

        // Mapping with patch mapping where possible

            //- Map field from default region to all non-default regions
            void map() const
            {
                forAll(mesh().regionNames(), regionI)
                {
                    word regionNameI = mesh().regionName(regionI);

                    if (regionNameI != polyMesh::defaultRegion)
                    {
                        map(regionI);
                    }
                }
            };

            //- Map field from default region to non-default region
            void map(label regionI) const
            {
                mapInteralField(regionI);
                mapBoundaryField(regionI);
            };

            //- Map field from default region to non-default region
            void map(const word& regionName) const
            {
                map(regions()[regionName]);
            };

            //- Map field from all non-default regions to default region
            void rmap() const
            {
                forAll(mesh().regionNames(), regionI)
                {
                    word regionNameI = mesh().regionName(regionI);

                    if (regionNameI != polyMesh::defaultRegion)
                    {
                        rmap(regionI);
                    }
                }
            };

            //- Map field from non-default region to default region
            void rmap(label regionI) const
            {
                rmapInteralField(regionI);
            };

            //- Map field from non-default region to default region
            void rmap(const word& regionName) const
            {
                rmap(regions()[regionName]);
            };

        // Mapping with patch internal values

            //- Map/Interpolate field from default region to all non-default regions
            void mapCopyInternal() const
            {
                forAll(mesh().regionNames(), regionI)
                {
                    word regionNameI = mesh().regionName(regionI);

                    if (regionNameI != polyMesh::defaultRegion)
                    {
                        mapCopyInternal(regionI);
                    }
                }
            };

            //- Map/Interpolate field from default region to non-default region
            void mapCopyInternal(label regionI) const
            {
                map(regionI);
                copyInternalBoundaryField(regionI);
            };

            //- Map/Interpolate field from default region to non-default region
            void mapCopyInternal(const word& regionName) const
            {
                mapCopyInternal(regions()[regionName]);
            };

        // Mapping with patch interpolation

            //- Map/Interpolate field from default region to all non-default regions
            void mapInterpolate() const
            {
                forAll(mesh().regionNames(), regionI)
                {
                    word regionNameI = mesh().regionName(regionI);

                    if (regionNameI != polyMesh::defaultRegion)
                    {
                        mapInterpolate(regionI);
                    }
                }
            };

            //- Map/Interpolate field from default region to non-default region
            void mapInterpolate(label regionI) const
            {
                map(regionI);
                interpolateBoundaryField(regionI);
            };

            //- Map/Interpolate field from default region to non-default region
            void mapInterpolate(const word& regionName) const
            {
                mapInterpolate(regions()[regionName]);
            };

        // Mapping with patch extrapolation

            //- Map/Extrapolate field from default region to all non-default regions
            void mapExtrapolate() const
            {
                forAll(mesh().regionNames(), regionI)
                {
                    word regionNameI = mesh().regionName(regionI);

                    if (regionNameI != polyMesh::defaultRegion)
                    {
                        mapExtrapolate(regionI);
                    }
                }
            };

            //- Map/Extrapolate field from default region to non-default region
            void mapExtrapolate(label regionI) const
            {
                map(regionI);
                extrapolateBoundaryField(regionI);
            };

            //- Map/Extrapolate field from default region to non-default region
            void mapExtrapolate(const word& regionName) const
            {
                mapExtrapolate(regions()[regionName]);
            };

    // Member Operators

        //- Field access operator by index
        inline GeometricField<Type, PatchField, GeoMesh>& operator[]
        (
            label regionI
        ) const
        {
            return fieldPtrs_[regionI];
        }

        //- Field access operator by regionName
        inline GeometricField<Type, PatchField, GeoMesh>& operator[]
        (
            const word& regionName
        ) const
        {
            return fieldPtrs_[regions()[regionName]];
        }

        void operator=
        (
            const regionGeometricField
            <
                Type, PatchField, GeoMesh, RegionGeoMesh
            >&
        );

        void operator=
        (
            const tmp
            <
                regionGeometricField
                <
                    Type, PatchField, GeoMesh, RegionGeoMesh
                >
            >&
        );

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "regionGeometricField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
