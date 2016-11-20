
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
//         volVectorField& ARe = globalStorage.ARe()[CONDUCTOR];
//         volVectorField& AIm = globalStorage.AIm()[CONDUCTOR];

// TODO: Use relaxation?
//         ARe.storePrevIter();
//         AIm.storePrevIter();

        // Update A
//         globalStorage.ARe().mapExtrapolate(CONDUCTOR);
//         globalStorage.AIm().mapExtrapolate(CONDUCTOR);
        globalStorage.ARe().mapInterpolate(CONDUCTOR);
        globalStorage.AIm().mapInterpolate(CONDUCTOR);

// TODO: Use relaxation?
//         ARe.relax(0.5);
//         AIm.relax(0.5);

        // Update rate of change of flux of A
        if (control.meshIs3D())
        {
            storage.ddtARePhi() =  omega0*(fvc::interpolate(globalStorage.ARe()[CONDUCTOR]) & mesh.Sf());
            storage.ddtAImPhi() = -omega0*(fvc::interpolate(globalStorage.AIm()[CONDUCTOR]) & mesh.Sf());

            if (globalSettings.biotSavart)
            {
                storage.ddtARePhi() += omega0*(fvc::interpolate(globalStorage.A0Re()[CONDUCTOR]) & mesh.Sf());
                storage.ddtAImPhi() -= omega0*(fvc::interpolate(globalStorage.A0Im()[CONDUCTOR]) & mesh.Sf());
            }
        }
    }
