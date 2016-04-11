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
#include "polyBoundaryMesh.H"
#include "foamTime.H"
#include "wordIOList.H"

// TODO [High]: Finish parallel face mapping
// TODO [High]: faceRegionAddressing and faceMap does not seem to be the same (flip-sign!!!)
// TODO [High]: How to deal with faceRegionAddressing and faceMap, properly?

// TODO [High]: Sort that parallel clutter out (remove local maps
//              and return maps based on if)!!!

// TODO [High]: Remove heaps of duplicated code

// TODO [High]: Remove copies of makeRegionPolyMeshData functions

// TODO [High]: Allow default region to be absent. WE NEED A LOT OF CHECKS!!!

// TODO [High]: Check if all region names are unique

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
private:

    // Private data

        //- Meshes data
        mutable List<polyMesh*> polyMeshes_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionPolyMesh(const regionPolyMesh&);

        //- Disallow default bitwise assignment
        void operator=(const regionPolyMesh&);

        // Helper functions for constructor

            //- Create new mesh for region with index regionI
            polyMesh* newMesh(const label& regionI) const;


protected:

    // Protected data

        //- Reference to time
        const Time& time_;

        //- Processor mesh?
        bool parallel_;

        //- Global region count
        mutable label size_;

        //- List of all region names
        mutable wordList regionNames_;

        //- Processor mesh: Was the mesh split into regions
        //  AFTER parallel decomposition of processors?
        mutable bool parallelSplitRegions_;

        // - Is mesh initialized?
        mutable bool initialized_;

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

            //- Read region names from file
            wordList readRegionNames() const;

            //- Resize lists storing data
            virtual void resizeLists() const;

            //- Init regions
            virtual void initMeshes(const wordList& regionNames) const;

            //- Set parallel split type
            void setParallelSplitRegions() const;

        // Make demand-driven data

            //- Region maps
            void makeCellRegionMap(const label& regionI) const;
            void makePointRegionMap(const label& regionI) const;
            void makeFaceRegionMap(const label& regionI) const;

            //- Processor maps
            void makeCellProcMap(const label& regionI) const;
            void makePointProcMap(const label& regionI) const;
// TODO            void makeFaceProcMap(const label& regionI) const;

            //- Region-processor maps
            void makeCellRegionProcMap(const label& regionI) const;
            void makePointRegionProcMap(const label& regionI) const;
// TODO            void makeFaceRegionProcMap(const label& regionI) const;

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


public:

    // Public typedefs

        typedef polyMesh Mesh;
        typedef polyBoundaryMesh BoundaryMesh;

    //- Runtime type information
    TypeName("regionPolyMesh");

    // Constructors

        //- Construct from time
        regionPolyMesh
        (
            const Time& runTime,
            bool init = true
        );

        //- Construct from time and list of region names
        regionPolyMesh
        (
            const Time& runTime,
            const wordList& regionNames,
            bool init = true
        );

    // Destructor

        virtual ~regionPolyMesh() {};

    // Member Functions

        // Initialization

            //- Init meshes
            void init(const wordList& regionNames) const
            {
                if (!initialized()) initMeshes(regionNames);
            };

            inline bool initialized() const
            {
                return initialized_;
            }

        // Access

            //- Return reference to regionPolyMesh
            inline const regionPolyMesh& rpmesh() const
            {
                return *this;
            };

            //- Return const reference to time
            inline const Time& time() const
            {
                return time_;
            };

            //- Return const reference to global region count
            inline const label& size() const
            {
                return size_;
            };

            //- Return const reference to the list of all region names
            inline const List<word>& regionNames() const
            {
                return regionNames_;
            };

            //- Return name of region with index regionI
            inline word regionName(const label& regionI) const
            {
                return regionNames_[regionI];
            }

// TODO: Use hash table
            //- Return region index for region with regionName
            label regionIndex
            (
                const word& regionName = polyMesh::defaultRegion
            ) const;

            //- Return whether this is a processor mesh
            inline const bool& parallel() const
            {
                return parallel_;
            }

            //- Return whether mesh was split into regions
            //  AFTER parallel decomposition of processors
            inline const bool& parallelSplitRegions() const
            {
                return parallelSplitRegions_;
            }

            //- Return reference to mesh
            inline polyMesh& mesh(const label& regionI) const
            {
                return operator[](regionI);
            }

            //- Return reference to mesh
            inline polyMesh& mesh
            (
                const word& regionName = polyMesh::defaultRegion
            ) const
            {
                return operator[](regionIndex(regionName));
            }

            //- Return the object registry
            inline const objectRegistry& thisDb() const
            {
                return *this;
            }

        // Active mesh data

            //- Return the vector of geometric directions in mesh
            inline const Vector<label>& geometricD(const label& regionI) const
            {
                return operator[](regionI).geometricD();
            };

            //- Return the vector of geometric directions in mesh
            inline const Vector<label>& geometricD
            (
                const word& regionName = polyMesh::defaultRegion
            ) const
            {
                return operator[](regionIndex(regionName)).geometricD();
            };

            //- Return the number of valid geometric dimensions in the mesh
            inline label nGeometricD(const label& regionI) const
            {
                return operator[](regionI).nGeometricD();
            };

            //- Return the number of valid geometric dimensions in the mesh
            inline label nGeometricD
            (
                const word& regionName = polyMesh::defaultRegion
            ) const
            {
                return operator[](regionIndex(regionName)).nGeometricD();
            };

            //- Return the vector of solved-for directions in mesh
            inline const Vector<label>& solutionD(const label& regionI) const
            {
                return operator[](regionI).solutionD();
            };

            //- Return the vector of solved-for directions in mesh
            inline const Vector<label>& solutionD
            (
                const word& regionName = polyMesh::defaultRegion
            ) const
            {
                return operator[](regionIndex(regionName)).solutionD();
            };

            //- Return the number of valid solved-for dimensions in the mesh
            inline label nSolutionD(const label& regionI) const
            {
                return operator[](regionI).nSolutionD();
            };

            //- Return the number of valid solved-for dimensions in the mesh
            inline label nSolutionD
            (
                const word& regionName = polyMesh::defaultRegion
            ) const
            {
                return operator[](regionIndex(regionName)).nSolutionD();
            };

        // Mapping

            //- Map mesh points from default region and replace them
            //  in a given point field of non-default region
            void map
            (
                pointField& givenPoints,
                const label& regionI,
                const word& patchName = ""
            ) const;

            //- Map mesh points from default region and replace them
            //  in a given point field of non-default region
            void map
            (
                pointField& givenPoints,
                const word& regionName,
                const word& patchName = ""
            ) const;

            //- Map mesh points from default region
            //  to a new point field for non-default region
            tmp<pointField> map
            (
                const label& regionI,
                const word& patchName = ""
            ) const;

            //- Map mesh points from default region
            //  to a new point field for non-default region
            tmp<pointField> map
            (
                const word& regionName,
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

        //- Mesh access operator by index
        virtual polyMesh& operator[](const label& regionI) const
        {
            return *polyMeshes_[regionI];
        }

        //- Mesh access operator by regionName
        inline polyMesh& operator[](const word& regionName) const
        {
            return operator[](regionIndex(regionName));
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
