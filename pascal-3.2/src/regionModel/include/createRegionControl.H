//
// createRegionControl.H
// ~~~~~~~~~~~~~~~~~~

    Foam::Info << "Create region control"
        << Foam::nl << Foam::endl;

    Foam::regionControl
    <
        Foam::regionGeoMesh <typeof(mesh_)>
    > control_(mesh_);
