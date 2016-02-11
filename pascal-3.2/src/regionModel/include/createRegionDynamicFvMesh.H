//
// createRegionDynamicFvMesh.H
// ~~~~~~~~~~~~~~~~~~

    Foam::Info
        << "Create dynamic region mesh for time = "
        << runTime.timeName() << endl;

    Foam::regionDynamicFvMesh mesh_(runTime);

    Foam::Info << "Regions: ";

    forAll (mesh_.regionNames(), regionI)
    {
        Foam::Info << mesh_.regionName(regionI);
        if (regionI < mesh_.regionNames().size()-1)
        {
            Foam::Info << ", ";
        }
    }

    Foam::Info << Foam::nl << Foam::endl;
