
    // Region mesh
    Info << "Create region mesh";

    regionFvMesh mesh_(runTime);

    Info << " (";

    forAll (mesh_.regionNames(), regionI)
    {
        Info << mesh_.regionName(regionI);
        if (regionI < mesh_.regionNames().size()-1)
        {
            Info << ", ";
        }
    }

    Info << ")" << nl << endl;
