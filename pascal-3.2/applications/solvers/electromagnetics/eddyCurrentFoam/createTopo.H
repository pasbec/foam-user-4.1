
    // TODO: What if conductor is in contact with boundary?
    // TODO: What happens on parallel?

    // Create set for conducting/non-conducting cells (cn):
    labelHashSet cCellsHashSet;
    labelHashSet nCellsHashSet;

    // Create sets for conducting/non-conducting faces (cn):
    // cnPN) Owner is in conducting region and neighbour
    //       lies in non-conducting region
    // cnNP) Owner is in non-conducting region and
    //       neighbour lies in conducting region
    labelHashSet cnFacesHashSet;
    labelHashSet cnPNfacesHashSet;
    labelHashSet cnNPfacesHashSet;

    label refFace = -1;
    label cRefCell = -1;
    label nRefCell = -1;
    label csRefCell = -1;

    // Assemble topologies
    {
        forAll(mesh.cells(), cellI)
        {
            if (sigma.internalField()[cellI] > SMALL)
            {
                cCellsHashSet.insert(cellI);
            }
            else
            {
                nCellsHashSet.insert(cellI);
            }
        }

        // Mesh references
        const unallocLabelList& own = mesh.owner();
        const unallocLabelList& nei = mesh.neighbour();

        // Find and insert faces
        forAll (own, facei)
        {
            if
            (
                (sigma.internalField()[own[facei]] > SMALL)
             && (sigma.internalField()[nei[facei]] < SMALL)
            )
            {
                cnFacesHashSet.insert(facei);
                cnPNfacesHashSet.insert(facei);
            }

            if
            (
                (sigma.internalField()[own[facei]] < SMALL)
             && (sigma.internalField()[nei[facei]] > SMALL)
            )
            {
                cnFacesHashSet.insert(facei);
                cnNPfacesHashSet.insert(facei);
            }
        }

//         // Set reference settings
//         // FIXME: Refpoint problems!
//         refFace = cnPNfacesHashSet.toc()[0];
//         cRefCell = own[refFace];
//         nRefCell = nei[refFace];

        // Set reference settings
        // TEST
        // FIXME: Refpoint problems!
        cRefCell = mesh.findCell(point(0.00, 0.00, 0.01));
        nRefCell = mesh.findCell(point(0.12, 0.00, 0.01));
    }

    // Create table of content lists
    labelList cCells = cCellsHashSet.toc();
    labelList nCells = nCellsHashSet.toc();

    labelList cnFaces = cnFacesHashSet.toc();
    labelList cnPNfaces = cnPNfacesHashSet.toc();
    labelList cnNPfaces = cnNPfacesHashSet.toc();

    // Create topoSets
    cellSet cCellSet
    (
        mesh,
        "cellSet_solver_c",
        cCellsHashSet
    ); cCellSet.write();


    cellSet nCellSet
    (
        mesh,
        "cellSet_solver_n",
        nCellsHashSet
    ); nCellSet.write();

    faceSet cnFaceSet
    (
        mesh,
        "faceSet_solver_cn",
        cnFacesHashSet
    ); cnFaceSet.write();

    faceSet cnPNfaceSet
    (
        mesh,
        "faceSet_solver_cn_PN",
        cnPNfacesHashSet
    ); cnPNfaceSet.write();

    faceSet cnNPfaceSet
    (
        mesh,
        "faceSet_solver_cn_NP",
        cnNPfacesHashSet
    ); cnNPfaceSet.write();



    // Create submesh
    // FIXME: Only if subset-solver is active!
    fvMeshSubset cMeshSubSet
    (
        IOobject
        (
            "conductor",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        ),
        mesh
    );

    // Initialize large cell subset
    cMeshSubSet.setLargeCellSubset(cCellsHashSet);

    // Create sub-mesh
    fvMesh& cMesh = cMeshSubSet.subMesh();

    // Mappings
    labelList cMeshCellMap = cMeshSubSet.cellMap();
    labelList cMeshFaceMap = cMeshSubSet.faceMap();

    // Reverse map reference cell for conducting region
    forAll(cMesh.C(), celli)
    {
        if (cMeshCellMap[celli] == cRefCell)
        {
            csRefCell = celli;
        }
    }
