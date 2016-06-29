
    fileName cMeshPrefix = "carrierMesh";
    if (mesh.name() != polyMesh::defaultRegion)
    {
        cMeshPrefix = mesh.name() + "_" + cMeshPrefix;
    }

    fvMeshSubset cMeshSubset
    (
        IOobject
        (
            cMeshPrefix,
            mesh.pointsInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    cMeshSubset.setCellSubset(aBMdata.faCells);
    fvMesh& cMesh = cMeshSubset.subMesh();
    cMesh.write();
