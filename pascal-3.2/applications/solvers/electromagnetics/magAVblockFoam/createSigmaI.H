
    // Create pointer for modified sigma indicator
    volScalarField* sigmaIcPtr = NULL;
    volScalarField* sigmaInPtr = NULL;

    // Calculate sigma indicator
    {
        // Create field for sigma indicator
        sigmaIcPtr = new volScalarField
        (
            IOobject
            (
                "sigmaIc",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<scalar>(word(),dimless,pTraits<scalar>::zero),
            zeroGradientFvPatchField<scalar>::typeName
        );
        sigmaInPtr = new volScalarField
        (
            IOobject
            (
                "sigmaIn",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<scalar>(word(),dimless,pTraits<scalar>::zero),
            zeroGradientFvPatchField<scalar>::typeName
        );

        // Get references to sigmaI pointer
        volScalarField& sigmaIc = *sigmaIcPtr;
        volScalarField& sigmaIn = *sigmaInPtr;

        // Set sigma to one in conducting/non-conducting region
        forAll(cCells, celli)
        {
            label cellI = cCells[celli];
            sigmaIc[cellI] = pTraits<scalar>::one;
        }
        forAll(nCells, celli)
        {
            label cellI = nCells[celli];
            sigmaIn[cellI] = pTraits<scalar>::one;
        }

        // Write sigmaI
        sigmaIc.write();
        sigmaIn.write();
    }
