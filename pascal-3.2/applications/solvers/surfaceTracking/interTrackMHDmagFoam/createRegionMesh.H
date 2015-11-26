
    // Region IDs
    const label defaultRegionID = 0;
    const label fluidRegionID = 1;
    const label dynamicRegionID = 2;

    // Region name list
    wordList regionNames (3, "");
    regionNames[defaultRegionID] = polyMesh::defaultRegion;
    regionNames[fluidRegionID] = "region_fluid";
    regionNames[dynamicRegionID] = "region_dynamic";


    // Create region mesh
    Info << "Create region mesh" << endl << endl;;
    regionDynamicFvMesh mesh_(runTime, regionNames);

//     // Region point mesh
//     Info << "Create region point mesh"
//         << endl;
//     regionPointMesh pointMesh_(mesh_);
//     Info << endl;
