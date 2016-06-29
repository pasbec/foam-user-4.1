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

Application
    makeFaSubPolyMesh

Description
    A mesh generator for a polyMesh (without cells) from a triangulated faMesh

\*---------------------------------------------------------------------------*/

#include "objectRegistry.H"
#include "foamTime.H"
#include "argList.H"
#include "OSspecific.H"
#include "faMesh.H"
#include "fvMesh.H"
#include "cellSet.H"
#include "triSurface.H"
#include "fvMeshSubset.H"
#include "Xfer.H"

using namespace Foam;

struct meshData
{
    const fvMesh& mesh;
    int nInternalFaces;
    int nBoundaryFaces;
    int nActiveFaces;
    int nInactiveFaces;
    meshData
    (
        const fvMesh& mesh_
    )
    :
        mesh(mesh_),
        nInternalFaces(mesh.nInternalFaces()),
        nBoundaryFaces(-1),
        nActiveFaces(-1),
        nInactiveFaces(-1)
    {
        nBoundaryFaces = 0;
        forAll (mesh.boundaryMesh(), patchI)
        {
            nBoundaryFaces += mesh.boundaryMesh()[patchI].size();
        }

        nActiveFaces = nInternalFaces + nBoundaryFaces;

        nInactiveFaces = mesh.nFaces() - nActiveFaces;
    }
};

struct aMeshData
{
    const fvMesh& mesh;
    const faMesh& aMesh;
    const labelList& meshPoints;
    const labelList& faceLabels;
    labelList faceCells;
    cellSet faInternalCells;
    cellSet faBoundaryCells;
    cellSet faCells;
    aMeshData
    (
        const fvMesh& mesh_,
        const faMesh& aMesh_
    )
    :
        mesh(mesh_),
        aMesh(aMesh_),
        meshPoints(aMesh.patch().meshPoints()),
        faceLabels(aMesh.faceLabels()),
        faceCells(faceLabels.size(),-1),
        faInternalCells
        (
            mesh,
            "cellSet_faInternalCells",
            labelHashSet()
        ),
        faBoundaryCells
        (
            mesh,
            "cellSet_faBoundaryCells",
            labelHashSet()
        ),
        faCells
        (
            mesh,
            "cellSet_faCells",
            labelHashSet()
        )
    {
        // Calculate face cell addressing
        forAll (faceLabels, faceI)
        {
            label globalFaceI = faceLabels[faceI];

            label patchID = mesh.boundaryMesh().whichPatch(globalFaceI);
            label patchStart = mesh.boundaryMesh()[patchID].start();
            label patchFaceI = globalFaceI - patchStart;

            label patchFaceCellI =
                mesh.boundaryMesh()[patchID].faceCells()[patchFaceI];

            faceCells[faceI] = patchFaceCellI;
        }

        // Fill internal cell set
        faInternalCells.insert(faceCells);

        // Add all cells in contact with edges of the boundary
        // to the cell set faBoundaryCells
        forAll (aMesh.boundary(), patchI)
        {
            const faPatch& aPatch = aMesh.boundary()[patchI];

            // Get neighbour patch
            label ngbPolyPatchID = aPatch.ngbPolyPatchIndex();
            if (ngbPolyPatchID != -1)
            {
                labelList ngbPolyPatchFaces = aPatch.ngbPolyPatchFaces();

                forAll (aPatch, edgeI)
                {
                    label ngbGlobalFaceI = ngbPolyPatchFaces[edgeI];

                    label ngbPatchStart = mesh.boundaryMesh()[ngbPolyPatchID].start();
                    label ngbPatchFaceI = ngbGlobalFaceI - ngbPatchStart;

                    // ngbPatch cell
                    label ngpPatchFaceCellI =
                        mesh.boundaryMesh()[ngbPolyPatchID].faceCells()[ngbPatchFaceI];

                    faBoundaryCells.insert(ngpPatchFaceCellI);
                }
            }
        }

        // Remove all cells in contact with internal faces
        // from the cell set faBoundaryCells
        forAllConstIter (labelHashSet, faInternalCells, iter)
        {
            faBoundaryCells.unset(iter.key());
        }

        // Fill cell set
        faCells.insert(faInternalCells.toc());
        faCells.insert(faBoundaryCells.toc());

        // Write cellSets
        faInternalCells.Foam::regIOobject::write();
        faBoundaryCells.Foam::regIOobject::write();
        faCells.Foam::regIOobject::write();
    }
};

class aMeshTriangulation
{
private:
    const faMesh& aMesh_;
    pointField* points_;
    triFaceList* faces_;
    labelList* faceMap_;
    List<labelHashSet>* faceRmap_;
private:
    void calcPoints()
    {
        int nBasePoints = aMesh_.nPoints();
        int nCentrePoints = aMesh_.nFaces();
        int nPoints = nBasePoints + nCentrePoints;

        points_ = new pointField(nPoints, vector::zero);
        pointField& points = *points_;

        const pointField& basePoints = aMesh_.points();
        const pointField& centrePoints = aMesh_.patch().faceCentres();

        // Base vertices
        forAll(basePoints, basePointI)
        {
            label pointI = basePointI;

            points[pointI] = basePoints[pointI];
        }

        // Face centre vertices
        forAll(centrePoints, centrePointI)
        {
            points[nBasePoints + centrePointI] = centrePoints[centrePointI];
        }
    }

    void calcTriangulation()
    {

        int nBaseEdges = aMesh_.nEdges();
        int nBaseInternalEdges = aMesh_.nInternalEdges();
        int nFaces = nBaseEdges + nBaseInternalEdges;

        faces_ = new triFaceList(nFaces);
        triFaceList& faces = *faces_;

        faceMap_ = new labelList(nFaces, 0);
        labelList& faceMap = *faceMap_;

        faceRmap_ = new List<labelHashSet>(aMesh_.nFaces(), labelHashSet());
        List<labelHashSet>& faceRmap = *faceRmap_;

        const edgeList& baseEdges = aMesh_.edges();

        // Internal triangulation based on lower/upper addressing
        int nBasePoints = aMesh_.nPoints();
        const unallocLabelList& baseOwn = aMesh_.owner();
        const unallocLabelList& baseNei = aMesh_.neighbour();

        forAll(baseOwn, baseEdgeI)
        {
            edge e = baseEdges[baseEdgeI];

            // Owner side triangle
            label ownP1 = e.start();
            label ownP2 = e.end();
            label ownP3 = nBasePoints + baseOwn[baseEdgeI];

            label ownFaceI = 2*baseEdgeI;

            faces[ownFaceI] = triFace(ownP1, ownP2, ownP3);
            faceMap[ownFaceI] = baseOwn[baseEdgeI];
            faceRmap[baseOwn[baseEdgeI]].insert(ownFaceI);

            // Neighbour side triangle
            label neiP1 = e.end();
            label neiP2 = e.start();
            label neiP3 = nBasePoints + baseNei[baseEdgeI];

            label neiFaceI = 2*baseEdgeI + 1;

            faces[neiFaceI] = triFace(neiP1, neiP2, neiP3);
            faceMap[neiFaceI] = baseNei[baseEdgeI];
            faceRmap[baseNei[baseEdgeI]].insert(neiFaceI);
        }

        // Patch boundary triangulation
        const faPatchList& basePatches = aMesh_.boundary();
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
                    faceMap[pFaceI] = basePatchEdgeFaces[basePatchEdgeI];
                    faceRmap[basePatchEdgeFaces[basePatchEdgeI]].insert(pFaceI);
                }
            }
        }

        // Empty boundary triangulation
        int nBaseNonEmptyEdges = nBaseInternalEdges + nBasePatchEdges;
        int nBaseEmptyEdges = nBaseEdges - nBaseNonEmptyEdges;

        if (nBaseEmptyEdges > 0)
        {
            labelList::subList baseEmptyOwn
            (
                aMesh_.edgeOwner(),
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
                faceMap[eFaceI] = baseEmptyOwn[baseEmptyEdgeI];
                faceRmap[baseEmptyOwn[baseEmptyEdgeI]].insert(eFaceI);
            }
        }
    }
public:
    aMeshTriangulation
    (
        const faMesh& aMesh
    )
    :
        aMesh_(aMesh),
        points_(NULL),
        faces_(NULL),
        faceMap_(NULL),
        faceRmap_(NULL)
    {
        calcPoints();
        calcTriangulation();
    }

    ~aMeshTriangulation()
    {
        delete faceRmap_;
        delete faceMap_;
        delete faces_;
        delete points_;
    }

    const pointField& points() const
    {
        return *points_;
    }

    const triFaceList& faces() const
    {
        return *faces_;
    }

    const labelList& faceMap() const
    {
        return *faceMap_;
    }

    const List<labelHashSet>& faceRmap() const
    {
        return *faceRmap_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
#   include "addRegionOption.H"

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createNamedMesh.H"
#   include "createFaMesh.H"

    // Base mesh data
    meshData bMdata(mesh);

Info << "DEBUG | bMdata.nInternalFaces = " << bMdata.nInternalFaces << endl;
Info << "DEBUG | bMdata.nBoundaryFaces = " << bMdata.nBoundaryFaces << endl;
Info << "DEBUG | bMdata.nActiveFaces = " << bMdata.nActiveFaces << endl;
Info << "DEBUG | bMdata.nInactiveFaces = " << bMdata.nInactiveFaces << endl;
Info << endl;

    // Area base-mesh data
    aMeshData aBMdata(mesh, aMesh);

Info << "DEBUG | aBMdata.internalCells.size() = " << aBMdata.faInternalCells.size() << endl;
Info << "DEBUG | aBMdata.boundaryCells.size() = " << aBMdata.faBoundaryCells.size() << endl;
Info << "DEBUG | aBMdata.cells.size() = " << aBMdata.faCells.size() << endl;
Info << endl;

    // Create carrier-mesh
#   include "createCarrierMesh.H"

    // Carrier-mesh data
    meshData cMdata(cMesh);

Info << "DEBUG | cMdata.nInternalFaces = " << cMdata.nInternalFaces << endl;
Info << "DEBUG | cMdata.nBoundaryFaces = " << cMdata.nBoundaryFaces << endl;
Info << "DEBUG | cMdata.nActiveFaces = " << cMdata.nActiveFaces << endl;
Info << "DEBUG | cMdata.nInactiveFaces = " << cMdata.nInactiveFaces << endl;
Info << endl;

    // Area base-mesh triangulation
    aMeshTriangulation aBMtri(aMesh);

Info << "DEBUG | aBMtri.points().size() = " << aBMtri.points().size() << endl;
Info << "DEBUG | aBMtri.faces().size() = " << aBMtri.faces().size() << endl;
// Info << "DEBUG | aBMtri.points() = " << aBMtri.points() << endl;
// Info << "DEBUG | aBMtri.faces() = " << aBMtri.faces() << endl;
Info << endl;

// Points

    labelHashSet usedCMpointLabels;

    // Collect all base-mesh point labels from points in carrier-mesh
    forAll(cMeshSubset.pointMap(), pointI)
    {
            usedCMpointLabels.insert(cMeshSubset.pointMap()[pointI]);
    }

    // Remove all base-mesh point labels from area base-mesh
    forAll(aBMdata.meshPoints, pointI)
    {
            usedCMpointLabels.unset(aBMdata.meshPoints[pointI]);
    }

    // Init points
    int nPoints = usedCMpointLabels.size() + aBMtri.points().size();
    pointField points(nPoints, vector::zero);

    // Add used points of carrier-mesh
// TODO: Point order?
    labelList usedCMpointLabelList = usedCMpointLabels.toc();
    forAll(usedCMpointLabelList, pointI)
    {
        points[pointI] = mesh.points()[usedCMpointLabelList[pointI]];
    }

    // Add all points of triangulated area base-mesh
    forAll (aBMtri.points(), pointI)
    {
        points[usedCMpointLabels.size() + pointI] = aBMtri.points()[pointI];
    }

Info << "DEBUG | usedCMpointLabels.size() = " << usedCMpointLabels.size() << endl;
Info << "DEBUG | points.size() = " << points.size() << endl;
// Info << "DEBUG | usedCMpointLabels.toc() = " << usedCMpointLabels.toc() << endl;
// Info << "DEBUG | points = " << points << endl;
Info << endl;

// Faces

    labelHashSet usedCMfaceLabels;

    // Collect all base-mesh face labels from faces in carrier-mesh
    forAll(cMeshSubset.faceMap(), faceI)
    {
            usedCMfaceLabels.insert(cMeshSubset.faceMap()[faceI]);
    }

    // Remove all base-mesh face labels from area base-mesh
    forAll(aBMdata.faceLabels, faceI)
    {
            usedCMfaceLabels.unset(aBMdata.faceLabels[faceI]);
    }

    // Init faces
    int nFaces = usedCMfaceLabels.size() + aBMtri.faces().size();
    faceList faces(nFaces);

    // Create reverse pointMap from base-mesh to carrier-mesh
// TODO: Redcuce memory? At least tmp<> should be used here!
    labelList cMpointRmap(mesh.points().size(),-1);
    forAll (usedCMpointLabelList, pointI)
    {
        cMpointRmap[usedCMpointLabelList[pointI]] = pointI;
    }

    // Create reverse pointMap from base-mesh to area base-mesh
// TODO: Redcuce memory? At least tmp<> should be used here!
    labelList aBMpointRmap(mesh.points().size(),-1);
    forAll (aBMdata.meshPoints, pointI)
    {
        aBMpointRmap[aBMdata.meshPoints[pointI]] = usedCMpointLabels.size() + pointI;
    }

    // Add used faces of carrier-mesh
    labelList usedCMfaceLabelList = usedCMfaceLabels.toc();
    forAll(usedCMfaceLabelList, faceI)
    {
        face curFace = mesh.faces()[usedCMfaceLabelList[faceI]];

        // Map face vertices
        labelList pointLabels(curFace.size(), -1);
        forAll(pointLabels, pointI)
        {
            // Point label is part of usedCMpointLabels (point from carrier-mesh)
            if (usedCMpointLabels.found(curFace[pointI]))
            {
                pointLabels[pointI] = cMpointRmap[curFace[pointI]];
            }
            // Point label is NOT part of usedCMpointLabels (point from area base-mesh)
            else
            {
                pointLabels[pointI] = aBMpointRmap[curFace[pointI]];
            }
        }

        faces[faceI] = face(pointLabels);
    }

    // Add all faces of triangulated area base-mesh
    forAll (aBMtri.faces(), faceI)
    {
        face curFace = aBMtri.faces()[faceI];

        // Map face vertices
        labelList pointLabels(curFace.size(), -1);
        forAll(pointLabels, pointI)
        {
            pointLabels[pointI] = usedCMpointLabels.size() + curFace[pointI];
        }

        faces[usedCMfaceLabels.size() + faceI] = face(pointLabels);
    }

Info << "DEBUG | usedCMfaceLabels.size() = " << usedCMfaceLabels.size() << endl;
// Info << "DEBUG | usedCMfaceLabels.toc() = " << usedCMfaceLabels.toc() << endl;
Info << "DEBUG | faces.size() = " << faces.size() << endl;
// Info << "DEBUG | faces = " << faces << endl;
Info << endl;

labelHashSet allPointLabels;
forAll(faces, faceI)
{
    face curFace = faces[faceI];
    forAll(curFace, pointI)
    {
        allPointLabels.insert(curFace[pointI]);
    }
}
Info << "DEBUG | allPointLabels.size() = " << allPointLabels.size() << endl;
// Info << "DEBUG | allPointLabels.sortedToc() = " << allPointLabels.sortedToc() << endl;
Info << endl;

// TODO FIXME: Face ordering! (Internal/Boundary/Inactive)

// TODO: Cells

// TODO: Patches

        // Create polyMesh
        polyMesh sMesh
        (
            IOobject
            (
                mesh.name() + "_faSubMesh",
                mesh.pointsInstance(),
                runTime,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            xferCopy(points),
            xferCopy(faces),
            xferCopy(labelList(0)),
            xferCopy(labelList(0)),
            false
        );

// TODO: Use patches from cMesh as they may be different
//     // Create patch for area sub-mesh
//     polyPatch* aSubPatch = polyPatch::New
//         (
//             "patch",                  // type
//             "faMesh",                 // name
//             aBMtri.faces().size(),    // size
//             usedCMpointLabels.size(), // start
//             0,                        // index
//             sMesh.boundaryMesh()      // boundary mesh
//         ).ptr();
// 
//     // Add patches
//     List<polyPatch*> sPatches(1, aSubPatch);
// 
//     sMesh.addPatches(sPatches, true);

    // Write subMesh
    Info << "Write subMesh ... ";
    sMesh.write();
    Info << "Done" << endl;

//     // Write addressings
//     //
// 
//        Info << "Write faceSubAddressing ... ";
//         labelIOList ioFaceSubAddressing
//         (
//             IOobject
//             (
//                 "faceSubAddressing",
//                 mesh.facesInstance(),
//                 faMesh::meshSubDir,
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE
//             ),
//             faceSubAddressing
//         );
//         ioFaceSubAddressing.write();
//         Info << "Done" << endl;
// 
//         Info << "Write faceBaseAddressing ... ";
//         labelIOList ioFaceBaseAddressing
//         (
//             IOobject
//             (
//                 "faceBaseAddressing",
//                 subMesh.facesInstance(),
//                 faMesh::meshSubDir,
//                 subMesh,
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE
//             ),
//             faceBaseAddressing
//         );
//         ioFaceBaseAddressing.write();
//         Info << "Done" << endl;

    return(0);
}

// ************************************************************************* //

