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
    regionPolyMesh

Description
    ...

SourceFiles
    regionPolyMesh.C
    regionPolyMeshMapping.C
    makeRegionPolyMeshData.C

\*---------------------------------------------------------------------------*/

#ifndef regionPolyMesh_H
#define regionPolyMesh_H

#include "objectRegistry.H"
#include "polyMesh.H"
#include "foamTime.H"
#include "wordList.H"

// TODO [High]: Finish implementation of face mapping in parallel!

// TODO [High]: Sort that parallel clutter out (remove local maps and return based on if)!!!

// TODO [High]: Remove heaps of duplicated code

// TODO [High]: Remove copies of mappings

// TODO [High]: Check if all region names are unique

// TODO [High]: Implement forward (base->regionI) point mapping
//              ala map()-functions


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class regionPolyMesh Declaration
\*---------------------------------------------------------------------------*/

class regionPolyMesh
:
    public objectRegistry
{
protected:

    // Protected data

        //- Reference to time
        const Time& time_;

        //- Global region count
        label size_;

        //- List of all region names
        wordList regionNames_;

        //- Processor mesh?
        bool parallel_;

        //- Processor mesh: Was the mesh split into regions
        //  AFTER parallel decomposition of processors?
        bool parallelSplitRegions_;

        //- Meshes data
        List<polyMesh*> meshesData_;

        //- Meshes access pointers
        mutable List<polyMesh*> meshes_;

    // Protected demand-driven data

        //- Region cell mapping
        mutable List<labelIOList*> cellRegionMap_;

        //- Region point mapping
        mutable List<labelIOList*> pointRegionMap_;

        //- Region face mapping
        mutable List<labelIOList*> faceRegionMap_;

        //- Processor cell mapping
        mutable List<labelIOList*> cellProcMap_;

        //- Processor point mapping
        mutable List<labelIOList*> pointProcMap_;

        //- Processor face mapping
        mutable List<labelIOList*> faceProcMap_;

        //- Region-processor cell mapping
        //  (if split in order region->proc)
        mutable List<labelIOList*> cellRegionProcMap_;

        //- Region-processor point mapping
        //  (if split in order region->proc)
        mutable List<labelIOList*> pointRegionProcMap_;

        //- Region-processor face mapping
        //  (if split in order region->proc)
        mutable List<labelIOList*> faceRegionProcMap_;

        //- Local cell mapping
        mutable List<labelIOList*> cellMap_;

        //- Local point mapping
        mutable List<labelIOList*> pointMap_;

        //- Local face mapping
        mutable List<labelIOList*> faceMap_;

    // Protected Member Functions

        // Helper functions for constructor

            //- Set parallel split type
            void setParallelSplitRegions();

        // Make demand-driven data

            //- Region maps
            void makeCellRegionMap(const label& regionI) const;
            void makePointRegionMap(const label& regionI) const;
            void makeFaceRegionMap(const label& regionI) const;

            //- Processor maps
            void makeCellProcMap(const label& regionI) const;
            void makePointProcMap(const label& regionI) const;
//             void makeFaceProcMap(const label& regionI) const;

            //- Region-processor maps
            void makeCellRegionProcMap(const label& regionI) const;
            void makePointRegionProcMap(const label& regionI) const;
//             void makeFaceRegionProcMap(const label& regionI) const;

            //- Maps
            void makeCellMap(const label& regionI) const;
            void makePointMap(const label& regionI) const;
            void makeFaceMap(const label& regionI) const;

        // Helper functions for demand-driven data

            //- Create new region map
            //  for mesh of region with index regionI
            labelIOList* newRegionMap
            (
                const label& regionI,
                const word& topoType
            ) const;

            //- Create new processor map
            //  for mesh of region with index regionI
            labelIOList* newProcMap
            (
                const label& regionI,
                const word& topoType
            ) const;

            //- Create new region-processor map
            //  for mesh of region with index regionI
            labelIOList* newCellRegionProcMap
            (
                const label& regionI
            ) const;
            labelIOList* newPointRegionProcMap
            (
                const label& regionI
            ) const;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionPolyMesh(const regionPolyMesh&);

        //- Disallow default bitwise assignment
        void operator=(const regionPolyMesh&);

        // Helper functions for constructor

            //- Create new mesh for region with index regionI
            polyMesh* newMesh(const label& regionI) const;

public:

    //- Runtime type information
    TypeName("regionPolyMesh");

    //- Return the region data sub-directory name (usually "regionPolyMesh")
    static word regionPolyMeshSubDir;

    // Public typedefs

    // Constructors

        //- Construct from time and list of region names
        regionPolyMesh
        (
            const Time& runTime,
            const wordList& regionNames,
            bool init = true
        );

    // Destructor

        virtual ~regionPolyMesh();

    // Member Functions

        // Access

            //- Return reference to regionPolyMesh
            const regionPolyMesh& rpmesh() const
            {
                return *this;
            };

            //- Return const reference to time
            const Time& time() const
            {
                return time_;
            };

            //- Return const reference to global region count
            const label& size() const
            {
                return size_;
            };

            //- Return const reference to the list of all region names
            const List<word>& regionNames() const
            {
                return regionNames_;
            };

            //- Return name of region with index regionI
            word regionName(const label& regionI) const
            {
                return regionNames_[regionI];
            }

            //- Return whether this is a processor mesh
            const bool& parallel() const
            {
                return parallel_;
            }

            //- Return whether mesh was split into regions
            //  AFTER parallel decomposition of processors
            const bool& parallelSplitRegions() const
            {
                return parallelSplitRegions_;
            }

            //- Return region index for region with regionName
            label regionIndex(const word& regionName) const;

            //- Return reference to mesh access ist
            const List<polyMesh*>& meshes() const
            {
                return meshes_;
            };

            //- Return reference to mesh of default region
            const polyMesh& mesh() const
            {
                return *meshes_[regionIndex(polyMesh::defaultRegion)];
            }

            //- Return reference to mesh
            const polyMesh& mesh(const label& regionI) const
            {
                return *meshes_[regionI];
            }

            //- Return reference to mesh
            const polyMesh& mesh(const word& regionName) const
            {
                return *meshes_[regionIndex(regionName)];
            }

            //- Return the object registry
            const objectRegistry& thisDb() const
            {
                return *this;
            }

        // Mapping

            //- Map mesh points from all non-default region and replace them
            //  in a given point field of default region
            void rmap
            (
                pointField& givenPoints,
                const word& patchName = ""
            ) const;

            //- Map mesh points from non-default region and replace them
            //  in a given point field of default region
            void rmap
            (
                pointField& givenPoints,
                const label& regionI,
                const word& patchName = ""
            ) const;

            //- Map mesh points from non-default region and replace them
            //  in a given point field of default region
            void rmap
            (
                pointField& givenPoints,
                const word& regionName,
                const word& patchName = ""
            ) const;

            //- Map mesh points from all non-default regions
            //  to a new point field for default region
            tmp<pointField> rmap
            (
                const word& patchName = ""
            ) const;

            //- Map mesh points from non-default region
            //  to a new point field for default region
            tmp<pointField> rmap
            (
                const label& regionI,
                const word& patchName = ""
            ) const;

            //- Map mesh points from non-default region
            //  to a new point field for default region
            tmp<pointField> rmap
            (
                const word& regionName,
                const word& patchName = ""
            ) const;

            //- Return patch mapping between two neighbour-regions
            //  based on directMappedPatch
            labelListList patchMapDirectMapped
            (
                const label& fromRegionI,
                const label& toRegionI
            ) const;

            //- Return patch mapping between two neighbour-regions
            //  based on directMappedPatch
            labelListList patchMapDirectMapped
            (
                const word& fromRegionName,
                const word& toRegionName
            ) const;

        // Demand-driven access

            //- Return const reference to mesh cell mapping
            const labelIOList& cellMap(const label& regionI) const;

            //- Return const reference to mesh cell mapping
            const labelIOList& cellMap(const word& regionName) const;

            //- Return const reference to mesh point mapping
            const labelIOList& pointMap(const label& regionI) const;

            //- Return const reference to mesh point mapping
            const labelIOList& pointMap(const word& regionName) const;

            //- Return const reference to mesh face mapping
            const labelIOList& faceMap(const label& regionI) const;

            //- Return const reference to mesh face mapping
            const labelIOList& faceMap(const word& regionName) const;

    // Member operators

        //- Mesh access operator for default region (lookup)
        const polyMesh& operator()() const
        {
            return *meshes_[regionIndex(polyMesh::defaultRegion)];
        }

        //- Mesh access operator by index (direct)
        const polyMesh& operator[](const label& regionI) const
        {
            return *meshes_[regionI];
        }

        //- Mesh access operator by regionName (lookup)
        const polyMesh& operator()(const word& regionName) const
        {
            return *meshes_[regionIndex(regionName)];
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
