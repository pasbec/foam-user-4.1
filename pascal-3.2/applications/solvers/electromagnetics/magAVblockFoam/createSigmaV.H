
    // Do we need the electric scalar potential?
    bool requireV = false;

    // TODO: Sigma V only if sigma is NOT constant
    if (cCells.size() > 0)
    {
        requireV = true;
    }

    // Create pointer for modified sigma
    volScalarField* sigmaVPtr = NULL;

    // Calculate modified sigma if required
    if (requireV)
    {
        // Copy sigma and change bc to zero-gradient
        sigmaVPtr = new volScalarField
        (
            IOobject
            (
                "sigmaV",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            sigma,
            zeroGradientFvPatchField<scalar>::typeName
        );

        // Write sigmaV
        sigmaVPtr->write();

        // Get reference to sigmaV pointer
        volScalarField& sigmaV = *sigmaVPtr;

        // Assemble equation for modified sigma
        // Note: Even though sigmaV has zero-gradient
        //       bc, we will NOT set a reference as we
        //       will copy all real sigma values from
        //       the conducting region, below
        fvScalarMatrix sigmaVEqn(fvm::laplacian(sigmaV));

        // Create scalar field with sigma values of all
        // conducting cells
        scalarField cCellSigmaV(cCells.size(), 0.0);
        forAll(cCells, celli)
        {
            label cellI = cCells[celli];
            cCellSigmaV[celli] = sigmaV[cellI];
        }

        // Set values in sigmaVEqn for conducting cells
        // to real sigma values from the scalar field above
        sigmaVEqn.setValues(cCells, cCellSigmaV);

        // Solve for modified sigma
        sigmaVEqn.solve();
    }
