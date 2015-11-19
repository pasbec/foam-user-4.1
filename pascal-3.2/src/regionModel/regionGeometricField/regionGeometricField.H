
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

Class
    regionGeometricField

Description
    ...

SourceFiles
    regionGeometricField.C
    regionGeometricFieldMapping.C

\*---------------------------------------------------------------------------*/

#ifndef regionGeometricField_H
#define regionGeometricField_H

#include "regIOobject.H"
#include "GeometricFields.H"
#include "dimensionedType.H"

#include "regionGeoMesh.H"

// TODO [High]: Mapping should be based on virtual functions of field type cell, point, face, ...

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

    // Private type defs

    // Private data

        //- Global region count
        const label& size_;

        //- List of all region names
        const List<word>& regionNames_;

        //- Reference to regionMesh
        const RegionMesh& regionMesh_;

        //- Fields
        List<GeometricField<Type, PatchField, GeoMesh>*> fields_;

    // Private Member Functions

        // Helper functions for constructor

            //- Return IOobject for region with index regionI
            const IOobject
            regionIO
            (
                const label& regionI,
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
            const IOobject&,
            const RegionMesh&,
            const dimensioned<Type>&,
            const word& patchFieldType=PatchField<Type>::calculatedType()
        );

        //- Constructor given IOobject, mesh, dimensioned<Type>
        //  and patch types.
        regionGeometricField
        (
            const IOobject&,
            const RegionMesh&,
            const dimensioned<Type>&,
            const wordList& patchFieldTypes
        );

        //- Construct and read given IOobject
        regionGeometricField
        (
            const IOobject&,
            const RegionMesh&
        );

        //- Construct as copy resetting IO parameters
        regionGeometricField
        (
            const IOobject&,
            const regionGeometricField
            <
                Type, PatchField, GeoMesh,
                RegionGeoMesh
            >&
        );

    // Destructor

        virtual ~regionGeometricField();

    // Member Functions

        // Access

            //- Return const reference to region mesh
            const RegionMesh& mesh() const
            {
                return regionMesh_;
            };

            //- Return dimensions
            const dimensionSet& dimensions() const
            {
                return fields_
                [
                    mesh().regionIndex(polyMesh::defaultRegion)
                ]->dimensions();
            }

            //- Return reference to field of default region
            GeometricField<Type, PatchField, GeoMesh>&
            field() const
            {
                return *fields_
                [
                    mesh().regionIndex(polyMesh::defaultRegion)
                ];
            }

            //- Return reference to field
            GeometricField<Type, PatchField, GeoMesh>&
            field(const label& regionI) const
            {
                return *fields_[regionI];
            }

            //- Return reference to field
            GeometricField<Type, PatchField, GeoMesh>&
            field(const word& regionName) const
            {
                return *fields_[mesh().regionIndex(regionName)];
            }

            //- WriteData member function required by regIOobject
            bool writeData(Ostream&) const
            {
                return false;
            }

        // Mapping

            //- Map field from default region to all non-default regions
            void map() const;

            //- Map field from default region to non-default region
            void map(const label& regionI) const;

            //- Map field from default region to non-default region
            void map(const word& regionName) const;

            //- Map field from all non-default regions to default region
            void rmap() const;

            //- Map field from non-default region to default region
            void rmap(const label& regionI) const;

            //- Map field from non-default region to default region
            void rmap(const word& regionName) const;

    // Member operators

        //- Field access operator for default region (lookup)
        GeometricField<Type, PatchField, GeoMesh>&
        operator()() const
        {
            return field(polyMesh::defaultRegion);
        }

        //- Field access operator by index (direct)
        GeometricField<Type, PatchField, GeoMesh>&
        operator[](const label& regionI) const
        {
            return field(regionI);
        }

        //- Field access operator by regionName (lookup)
        GeometricField<Type, PatchField, GeoMesh>&
        operator()(const word& regionName) const
        {
            return field(regionName);
        }


        void operator=
        (
            const regionGeometricField
            <
                Type, PatchField, GeoMesh,
                RegionGeoMesh
            >&
        );
        void operator=
        (
            const tmp
            <
                regionGeometricField
                <
                    Type, PatchField, GeoMesh,
                    RegionGeoMesh
                >
            >&
        );

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
