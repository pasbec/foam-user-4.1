
    // Create pointer for modified sigma indicator
    volScalarField* indicatorVcPtr = NULL;
    volScalarField* indicatorVnPtr = NULL;

    // Calculate sigma indicator
    {
        // Create field for sigma indicator
        indicatorVcPtr = new volScalarField
        (
            IOobject
            (
                "indicatorVc",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<scalar>(word(),dimless,pTraits<scalar>::zero),
            zeroGradientFvPatchField<scalar>::typeName
        );
        indicatorVnPtr = new volScalarField
        (
            IOobject
            (
                "indicatorVn",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<scalar>(word(),dimless,pTraits<scalar>::zero),
            zeroGradientFvPatchField<scalar>::typeName
        );

        // Get references to indicator pointers
        volScalarField& indicatorVc = *indicatorVcPtr;
        volScalarField& indicatorVn = *indicatorVnPtr;

        // Set indicator to one in its corresponding
        // conducting (c) or non-conducting (n) region
        forAll(cCells, celli)
        {
            label cellI = cCells[celli];
            indicatorVc[cellI] = pTraits<scalar>::one;
        }
        forAll(nCells, celli)
        {
            label cellI = nCells[celli];
            indicatorVn[cellI] = pTraits<scalar>::one;
        }

        // Write sigmaI
        indicatorVc.write();
        indicatorVn.write();
    }
