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

// TODO [High]: Correct boundary conditions after construction in order to get
//              coupled patch values correct. This may be a problem!

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

        typedef GeometricField<Type, PatchField, GeoMesh> GeometricFieldType;

        typedef regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>
            RegionGeometricFieldType;


private:

    // Private data

        //- Reference to regionMesh
        const RegionMesh& regionMesh_;

        //- Field pointers
        mutable List<GeometricFieldType*> fieldPtrs_;

        //- Field active state
        mutable List<bool> fieldActive_;

        //- Field link state
        mutable List<bool> fieldLinked_;

    // Private Member Functions

        // Helper functions for constructor

            //- Return IOobject for region with index regionI
            const IOobject regionIOobject
            (
                label regionI,
                const IOobject& IOo
            ) const;

            //- Try to link field pointer
            bool linkFieldPtr
            (
                label regionI,
                const IOobject& IOo
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
        regionGeometricField
        (
            const IOobject& IOo,
            const RegionMesh& regionMesh,
            const dimensioned<Type>& dim,
            const HashTable<IOobject> IOoOverride = HashTable<IOobject>()
        );

        //- Constructor given IOobject, mesh, dimensioned<Type>
        //  and patch type.
        regionGeometricField
        (
            const IOobject& IOo,
            const RegionMesh& regionMesh,
            const dimensioned<Type>& dim,
            const word& patchFieldType,
            const HashTable<IOobject> IOoOverride = HashTable<IOobject>()
        );

        //- Constructor given IOobject, mesh, dimensioned<Type>
        //  and patch types.
        regionGeometricField
        (
            const IOobject& IOo,
            const RegionMesh& regionMesh,
            const dimensioned<Type>& dim,
            const wordList& patchFieldTypes,
            const HashTable<IOobject> IOoOverride = HashTable<IOobject>()
        );

        //- Construct as copy resetting IO parameters
        regionGeometricField
        (
            const IOobject& IOo,
            const RegionGeometricFieldType& rgf,
            const HashTable<IOobject> IOoOverride = HashTable<IOobject>()
        );


    // Constructor link-wrapper

        //- Constructor given IOobject, mesh, dimensioned<Type>
        static regionGeometricField* LinkOrNew
        (
            const IOobject& IOo,
            const RegionMesh& regionMesh,
            const dimensioned<Type>& dim,
            const HashTable<IOobject> IOoOverride = HashTable<IOobject>()
        );

        //- Constructor given IOobject, mesh, dimensioned<Type>
        //  and patch type.
        static regionGeometricField* LinkOrNew
        (
            const IOobject& IOo,
            const RegionMesh& regionMesh,
            const dimensioned<Type>& dim,
            const word& patchFieldType,
            const HashTable<IOobject> IOoOverride = HashTable<IOobject>()
        );

        //- Constructor given IOobject, mesh, dimensioned<Type>
        //  and patch types.
        static regionGeometricField* LinkOrNew
        (
            const IOobject& IOo,
            const RegionMesh& regionMesh,
            const dimensioned<Type>& dim,
            const wordList& patchFieldTypes,
            const HashTable<IOobject> IOoOverride = HashTable<IOobject>()
        );

        //- Construct as copy resetting IO parameters
        static regionGeometricField* LinkOrNew
        (
            const IOobject& IOo,
            const RegionGeometricFieldType& rgf,
            const HashTable<IOobject> IOoOverride = HashTable<IOobject>()
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
            void mapBoundaryBoundaryField(label regionI) const
            ;
//             {
//                 notImplemented(type() + "::mapBoundaryBoundaryField()");
//             };

//- TODO: Boundary field interpolation
//- TODO: Parallel?
            void mapInterpolateBoundaryInternalField(label regionI) const
            ;
//             {
//                 notImplemented(type() + "::mapInterpolateBoundaryInternalField()");
//             };

            void mapExtrapolateBoundaryInternalField(label regionI) const
            {
                notImplemented(type() + "::mapExtrapolateBoundaryInternalField()");
            };

        // Mapping

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
                mapBoundaryBoundaryField(regionI);
                mapInterpolateBoundaryInternalField(regionI);
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
                mapBoundaryBoundaryField(regionI);
                mapExtrapolateBoundaryInternalField(regionI);
            };

            //- Map/Extrapolate field from default region to non-default region
            void mapExtrapolate(const word& regionName) const
            {
                mapExtrapolate(regions()[regionName]);
            };

    // Member Operators

        //- Field access operator by index
        inline GeometricFieldType& operator[]
        (
            label regionI
        ) const
        {
            return *fieldPtrs_[regionI];
        };

        //- Field access operator by regionName
        inline GeometricFieldType& operator[]
        (
            const word& regionName
        ) const
        {
            return *fieldPtrs_[regions()[regionName]];
        };

        //- Field access operator by index
        inline const GeometricFieldType* operator()
        (
            label regionI
        ) const
        {
            return fieldPtrs_[regionI];
        };

        //- Field access operator by regionName
        inline const GeometricFieldType* operator()
        (
            const word& regionName
        ) const
        {
            return fieldPtrs_[regions()[regionName]];
        };

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
