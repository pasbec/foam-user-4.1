
    // Update gradient of V in default region
    // from data in conductor region
    if (control.meshIs3D())
    {
        SM_REGIONSCOPE(CONDUCTOR);

        if (Control::debug)
        {
            Info<< Control::typeName << " | VGradUpdate.H : "
                << "Update gradient of V."
                << endl;
        }

        // Update gradient of V
        globalStorage.VReGrad().rmap(CONDUCTOR);
        globalStorage.VImGrad().rmap(CONDUCTOR);
    }
