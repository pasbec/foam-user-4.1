
    // Create pointer for modified sigma
    volScalarField* sigmaVPtr = NULL;

    // Calculate modified sigma if required
//     if (requireV)
    {
        // Create artificial electrical conductivity. This
        // is done in order to solve the V equation for the
        // whole domain. We will copy sigma and change bc
        // to zero-gradient
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

        // Get reference to sigmaV pointer
        volScalarField& sigmaV = *sigmaVPtr;

        // Assemble equation for modified sigma
        // Note: Even though sigmaV has zero-gradient
        //       bc, we will NOT set a reference as we
        //       will copy and set all real sigma values
        //       from the conducting region
        fvScalarMatrix sigmaVEqn
        (
            fvm::laplacian(sigmaV)
        );

        // Copy and set sigma values in sigmaVEqn for
        // conducting region
        scalarField cCellSigmaV(cCells.size(), 0.0);
        forAll(cCells, celli)
        {
            label cellI = cCells[celli];
            cCellSigmaV[celli] = sigmaV[cellI];
        }
        sigmaVEqn.setValues(cCells, cCellSigmaV);

        // Solve for modified sigma
        sigmaVEqn.solve();

        // Write sigmaV
        sigmaV.write();
    }
