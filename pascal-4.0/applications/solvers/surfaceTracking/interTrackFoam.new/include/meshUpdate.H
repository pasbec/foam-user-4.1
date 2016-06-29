
    {
        trackedSurface& interface = storage.interface();

        interface.updateMesh(); // abs phi
        interface.updateDisplacementDirections();

        interface.predictPoints(); // abs phi
    }
