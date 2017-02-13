
    // Update A in conductor region
    // from data in default region
    {
        SM_REGIONSCOPE(CONDUCTOR);

        if (Control::debug)
        {
            Info<< Control::typeName << " | AUpdate.H : "
                << "Update A in conductor region."
                << endl;
        }

// TODO: Use relaxation?
        globalStorage.ARe()[CONDUCTOR].storePrevIter();
        globalStorage.ARe()[CONDUCTOR].storePrevIter();

        // Update A
        globalStorage.ARe().mapInterpolate(CONDUCTOR);
        globalStorage.AIm().mapInterpolate(CONDUCTOR);

// TODO: Use relaxation?
        const scalar AReRf = globalMesh[DEFAULT].solutionDict().fieldRelaxationFactor(globalStorage.ARe()[DEFAULT].name());
        const scalar AImRf = globalMesh[DEFAULT].solutionDict().fieldRelaxationFactor(globalStorage.AIm()[DEFAULT].name());
        globalStorage.ARe()[CONDUCTOR].relax(AReRf);
        globalStorage.ARe()[CONDUCTOR].relax(AImRf);

        // Update rate of change of flux of A
        if (control.meshIs3D())
        {
// TEST
// TODO: Interpolation with rMur?
            storage.ddtARePhi() =  omega0*(fvc::interpolate(globalStorage.ARe()[CONDUCTOR]) & mesh.Sf());
            storage.ddtAImPhi() = -omega0*(fvc::interpolate(globalStorage.AIm()[CONDUCTOR]) & mesh.Sf());

            if (globalSettings.biotSavart)
            {
                storage.ddtARePhi() += omega0*(fvc::interpolate(globalStorage.A0Re()[CONDUCTOR]) & mesh.Sf());
                storage.ddtAImPhi() -= omega0*(fvc::interpolate(globalStorage.A0Im()[CONDUCTOR]) & mesh.Sf());
            }
        }
    }
