//
// createRegionMesh.H
// ~~~~~~~~~~~~~~~~~~

    Foam::Info
        << "Create region mesh for time = "
        << runTime.timeName();

    Foam::regionFvMesh mesh_(runTime);

    Foam::Info << " (";

    forAll (mesh_.regionNames(), regionI)
    {
        Foam::Info << mesh_.regionName(regionI);
        if (regionI < mesh_.regionNames().size()-1)
        {
            Foam::Info << ", ";
        }
    }

    Foam::Info << ")" << Foam::nl << Foam::endl;
