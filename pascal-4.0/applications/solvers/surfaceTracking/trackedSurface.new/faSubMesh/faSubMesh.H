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
    faSubMesh

Description
    Implementation of a triangulated sub-mesh of a given base finite-area
    mesh

SourceFiles
    faSubMesh.C
    faSubMeshDemandDrivenData.C

\*---------------------------------------------------------------------------*/

#ifndef faSubMesh_H
#define faSubMesh_H

#include "fvCFD.H"
#include "polyMesh.H"
#include "polyPatch.H"
#include "faCFD.H"
#include "faMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class faSubMesh Declaration
\*---------------------------------------------------------------------------*/

class faSubMesh
{
    // Private data

        //- Name
        const word name_;

        //- Base faMesh
        const polyMesh& basePolyMesh_;

        //- Base faMesh
        const faMesh& baseAreaMesh_;

        const pointField& splitFacePoints_;

    // Permanent data

        //- Sub polyPatch
        polyPatch* subPolyPatchPtr_;

        //- Sub polyMesh
        polyMesh* subPolyMeshPtr_;

        //- Sub faMesh
        faMesh* subAreaMeshPtr_;


    // Demand-driven data

        //- Face curvatures
        mutable areaScalarField* faceCurvaturesPtr_;


    // Static Private Data


    // Private Member Functions

        //- Disallow default bitwise copy construct
        faSubMesh(const faSubMesh&);

        //- Disallow default bitwise assignment
        void operator=(const faSubMesh&);

        // Helpers

            //- Make/Collect points for sub polyMesh
            tmp<pointField> calcSubPolyMeshPoints() const;

        //- Initialization

            //- Make sub polyMesh
            void initSubPolyMesh();

            //- Make sub faMesh
            void initSubAreaMesh();

        // Demand-driven data

            //- Make face curvatures
            void makeFaceCurvatures() const;

            //- Clear geomatry related demand-driven data
            void clearGeom() const;

            //- Clear all demand-driven data
            void clearOut() const;


public:

    // Declare name of the class and it's debug switch
    ClassName("faSubMesh");


    // Constructors

        faSubMesh
        (
            const word name,
            const faMesh& baseAreaMesh,
            const pointField& splitFacePoints
        );


    // Destructor

        ~faSubMesh();


    // Member Functions

        // Helpers


        // Database

        // Access

            //- Name
            const word name() const
            {
                return name_;
            }

            //- Return access to base polyMesh
            const polyMesh& basePolyMesh() const
            {
                return basePolyMesh_;
            }

            //- Return access to base faMesh
            const faMesh& baseAreaMesh() const
            {
                return baseAreaMesh_;
            }

            //- Return access to sub polyMesh
            polyMesh& subPolyMesh() const
            {
                return *subPolyMeshPtr_;
            }

            //- Return access to sub faMesh
            faMesh& subAreaMesh() const
            {
                return *subAreaMeshPtr_;
            }

            //- Return points to split faces of base mesh
            const pointField& splitFacePoints() const
            {
                return splitFacePoints_;
            }

        // Mesh motion and mophing

            //- Is sub polyMesh moving
            bool moving() const
            {
                return subPolyMesh().moving();
            }

            //- Move points and update faMesh
            bool movePoints() const
            {
                clearGeom();

                subPolyMesh().movePoints
                (
                    calcSubPolyMeshPoints()
                );

                return true;
            }

        // Mapping

        // Demand-driven data

            //- Return face curvatures calculated on sub faMesh
            //  interpolated to base faMesh
            const areaScalarField& faceCurvatures() const;

        //- Write mesh
//         virtual bool write() const;


    // Member Operators
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
