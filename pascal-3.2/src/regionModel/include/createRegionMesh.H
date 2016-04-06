//
// createRegionMesh.H
// ~~~~~~~~~~~~~~~~~~

    Foam::Info
        << "Create region mesh for time = "
        << runTime.timeName() << endl;

    Foam::regionFvMesh regionMesh(runTime);

    Foam::Info << "Regions: ";

    forAll (regionMesh.regionNames(), regionI)
    {
        Foam::Info << regionMesh.regionName(regionI);
        if (regionI < regionMesh.regionNames().size()-1)
        {
            Foam::Info << ", ";
        }
    }

    Foam::Info << Foam::nl << Foam::endl;
