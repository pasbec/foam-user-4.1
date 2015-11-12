
    // Region IDs
    const label defaultRegionID = 0;
    const label conductorRegionID = 1;


    // Region name list
    wordList regionNames (2, "");
    regionNames[defaultRegionID] = polyMesh::defaultRegion;
    regionNames[conductorRegionID] = "region_conductor";


    // Create region mesh
    Info << "Create region mesh"
        << endl;
    regionFvMesh mesh_(runTime, regionNames);
    Info << endl;
