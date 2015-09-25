
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
