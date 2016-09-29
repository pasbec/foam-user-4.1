
    {
        if (Control::debug)
        {
            Info << "interTrackApp::Control | meshUpdate.H : "
                << "Update mesh and predict interface points."
                << endl;
        }

        trackedSurface& interface = storage.interface();

        interface.updateMesh(); // abs phi

        interface.predictPoints(); // abs phi
    }
