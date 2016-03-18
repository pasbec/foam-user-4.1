/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description

\*---------------------------------------------------------------------------*/

#include "faSubMesh.H"

#include "triSurface.H"
#include "Xfer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(faSubMesh, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

tmp<pointField> faSubMesh::calcSubPolyMeshPoints() const
{

    int nBasePoints = baseAreaMesh().nPoints();
    int nBaseControlPoints = baseAreaMesh().nFaces();
    int nPoints = nBasePoints + nBaseControlPoints;

    const pointField& basePoints = baseAreaMesh().points();
    const pointField& baseControlPoints = splitFacePoints();

    tmp<pointField> tpoints
    (
        new pointField(nPoints)
    );
    pointField& points = tpoints();

    // Base vertices
    forAll(basePoints, basePointI)
    {
        label pointI = basePointI;

        points[pointI] = basePoints[pointI];
    }

    // Base control points
    forAll(baseControlPoints, baseControlPointI)
    {
        label pointI = nBasePoints + baseControlPointI;

        points[pointI] = baseControlPoints[baseControlPointI];
    }

    return tpoints;
}


void faSubMesh::initSubPolyMesh()
{
    // Collect triangulation points
    //

        pointField points = calcSubPolyMeshPoints();


    // Calc triangulation
    //

        int nBaseEdges = baseAreaMesh().nEdges();
        int nBaseInternalEdges = baseAreaMesh().nInternalEdges();

        int nFaces = nBaseEdges + nBaseInternalEdges;
        triFaceList faces(nFaces);

        const edgeList& baseEdges = baseAreaMesh().edges();


    // Internal triangulation
    // based on lower/upper addressing

        int nBasePoints = baseAreaMesh().nPoints();
        const unallocLabelList& baseOwn = baseAreaMesh().owner();
        const unallocLabelList& baseNei = baseAreaMesh().neighbour();

        forAll(baseOwn, baseEdgeI)
        {
            edge e = baseEdges[baseEdgeI];

            // Owner side triangle
            label ownP1 = e.start();
            label ownP2 = e.end();
            label ownP3 = nBasePoints + baseOwn[baseEdgeI];

            label ownFaceI = 2*baseEdgeI;

            faces[ownFaceI] = triFace(ownP1, ownP2, ownP3);

            // Neighbour side triangle
            label neiP1 = e.end();
            label neiP2 = e.start();
            label neiP3 = nBasePoints + baseNei[baseEdgeI];

            label neiFaceI = 2*baseEdgeI + 1;

            faces[neiFaceI] = triFace(neiP1, neiP2, neiP3);
        }


    // Patch boundary triangulation
    //

        const faPatchList& basePatches = baseAreaMesh().boundary();

        int nBasePatchEdges = 0;

        forAll(basePatches, basePatchI)
        {
            const faPatch& basePatch = basePatches[basePatchI];

            if (basePatch.size() > 0)
            {
                const unallocLabelList& basePatchEdgeFaces =
                    basePatch.edgeFaces();

                forAll(basePatch, basePatchEdgeI)
                {
                    // Start corresponding to current patch
                    label baseEdgeI = basePatch.start() + basePatchEdgeI;

                    edge e = baseEdges[baseEdgeI];

                    // Patch edge triangle
                    label pP1 = e.start();
                    label pP2 = e.end();
                    label pP3 =
                        nBasePoints + basePatchEdgeFaces[basePatchEdgeI];

                    label pFaceI =
                        2*nBaseInternalEdges + nBasePatchEdges++;

                    faces[pFaceI] = triFace(pP1, pP2, pP3);
                }
            }
        }


    // Empty boundary triangulation
    //

        int nBaseNonEmptyEdges = nBaseInternalEdges + nBasePatchEdges;
        int nBaseEmptyEdges = nBaseEdges - nBaseNonEmptyEdges;

        if (nBaseEmptyEdges > 0)
        {
            labelList::subList baseEmptyOwn
            (
                baseAreaMesh().edgeOwner(),
                nBaseEmptyEdges,
                nBaseNonEmptyEdges
            );

            forAll(baseEmptyOwn, baseEmptyEdgeI)
            {
                // Start after all non-empty edges
                label baseEdgeI = nBaseNonEmptyEdges + baseEmptyEdgeI;

                edge e = baseEdges[baseEdgeI];

                // Empty edge triangle
                label eP1 = e.start();
                label eP2 = e.end();
                label eP3 = nBasePoints + baseEmptyOwn[baseEmptyEdgeI];

                label eFaceI =
                    2*nBaseInternalEdges + nBasePatchEdges + baseEmptyEdgeI;

                faces[eFaceI] = triFace(eP1, eP2, eP3);
            }
        }


    // Create polyMesh
    //

        // Create mesh with points and faces
        faceList polyFaces(faces.size());

        forAll (polyFaces, faceI)
        {
            polyFaces[faceI] = faces[faceI];
        }

        subPolyMeshPtr_ = new polyMesh
        (
            IOobject
            (
                name(),
                basePolyMesh().pointsInstance(),
                basePolyMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            xferCopy(points),
            xferCopy(polyFaces),
            xferCopy(labelList(0)),
            xferCopy(labelList(0)),
            false
        );


    // Create/Add polyPatch
    //

        // Create one single patch containing all faces
        subPolyPatchPtr_ =
            polyPatch::New
            (
                "patch",                        // type
                "default",                      // name
                nFaces,                         // size
                0,                              // start
                0,                              // index
                subPolyMeshPtr_->boundaryMesh() // boundary mesh
            ).ptr();

        // Add patches
        List<polyPatch*> patches(1, subPolyPatchPtr_);

        subPolyMeshPtr_->addPatches(patches, true);
}


void faSubMesh::initSubAreaMesh()
{
    // Use first (and only) patch of sub polyMesh
    // to construct sub faMesh
    subAreaMeshPtr_ = new faMesh(subPolyMesh(), 0);
}


void faSubMesh::clearGeom() const
{
    deleteDemandDrivenData(faceCurvaturesPtr_);
}


void faSubMesh::clearOut() const
{
    clearGeom();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

faSubMesh::faSubMesh
(
    const word name,
    const faMesh& baseAreaMesh,
    const pointField& splitFacePoints
)
:
name_(name),
basePolyMesh_(baseAreaMesh.mesh()),
baseAreaMesh_(baseAreaMesh),
splitFacePoints_(splitFacePoints),
subPolyPatchPtr_(NULL),
subPolyMeshPtr_(NULL),
subAreaMeshPtr_(NULL),
faceCurvaturesPtr_(NULL)
{
    if (debug)
    {
        Info<< "faSubMesh::calcSubPolyMeshPoints() : "
            << "init polyMesh"
            << endl;
    }

    initSubPolyMesh();

    if (debug)
    {
        Info<< "faSubMesh::faSubMesh() : "
            << "init areaMesh"
            << endl;
    }

    initSubAreaMesh();
}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

faSubMesh::~faSubMesh()
{
    clearOut();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
