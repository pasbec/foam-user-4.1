
    // Create set for conducting cells
    labelHashSet cCellsHashSet;

    forAll(mesh.cells(), cellI)
    {
        if (sigma.internalField()[cellI] > SMALL)
        {
            cCellsHashSet.insert(cellI);
        }
    }

    labelList cCells = cCellsHashSet.toc();

    cellSet cCellSet
    (
        mesh,
        "cellSet_solver_c",
        cCellsHashSet
    ); cCellSet.write();



    // Create set for non-conducting cells
    labelHashSet nCellsHashSet;

    forAll(mesh.cells(), cellI)
    {
        if (sigma.internalField()[cellI] < SMALL)
        {
            nCellsHashSet.insert(cellI);
        }
    }

    labelList nCells = nCellsHashSet.toc();

    cellSet nCellSet
    (
        mesh,
        "cellSet_solver_n",
        nCellsHashSet
    ); nCellSet.write();



    // Create sets for conducting/non-conducting faces (cn):
    // cnPN) Owner is in conducting region and neighbour
    //       lies in non-conducting region
    // cnNP) Owner is in non-conducting region and
    //       neighbour lies in conducting region
    labelHashSet cnFacesHashSet;
    labelHashSet cnPNfacesHashSet;
    labelHashSet cnNPfacesHashSet;

    {
        // Mesh references
        const unallocLabelList& own = mesh.owner();
        const unallocLabelList& nei = mesh.neighbour();

        // Find and insert faces
        forAll (own, faceI)
        {
            if
            (
                (sigma.internalField()[own[faceI]] > SMALL)
             && (sigma.internalField()[nei[faceI]] < SMALL)
            )
            {
                cnFacesHashSet.insert(faceI);
                cnPNfacesHashSet.insert(faceI);
            }

            if
            (
                (sigma.internalField()[own[faceI]] < SMALL)
             && (sigma.internalField()[nei[faceI]] > SMALL)
            )
            {
                cnFacesHashSet.insert(faceI);
                cnNPfacesHashSet.insert(faceI);
            }
        }
    }

    labelList cnFaces = cnFacesHashSet.toc();

    faceSet cnFaceSet
    (
        mesh,
        "cellSet_solver_cn",
        cnFacesHashSet
    ); cnFaceSet.write();

    faceSet cnPNfaceSet
    (
        mesh,
        "cellSet_solver_cn_PN",
        cnPNfacesHashSet
    ); cnPNfaceSet.write();

    faceSet cnNPfaceSet
    (
        mesh,
        "cellSet_solver_cn_NP",
        cnNPfacesHashSet
    ); cnNPfaceSet.write();
