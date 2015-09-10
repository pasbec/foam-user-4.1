
    // Interpolate sigma to face centers
    surfaceScalarField sigmaf
    (
        "sigmaf",
        fvc::interpolate(sigma)
    );

    // Calculate gradient of alpha
    volVectorField gradAlpha
    (
        IOobject
        (
            "gradAlpha",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        omega * fvc::grad(sigma),
        fixedValueFvPatchVectorField::typeName
    );
    gradAlpha.boundaryField() == vector::zero;

    // Manipulate face interpolated sigma:
    //
    //  We will replace zero values with very small
    //  values at all faces in order to avoid a singular
    //  solution
    {
        // Parameters
        scalar sigmafNonConducting = SMALL;

        // Manipulate internal faces
        forAll(sigmaf.internalField(), faceI)
        {
            if
            (sigmaf.internalField()[faceI] < SMALL)
            {
                sigmaf.internalField()[faceI] = sigmafNonConducting;
            }
        }

        // Manipulate boundary faces
        sigmaf.boundaryField() == sigmafNonConducting;
    }

    // Reference settings
    label refCell = mesh.findCell(point(0.00, 0.00, 0.01));
    scalar refValue = 0.0;

    // Real part
    fvScalarMatrix VReEqn
    (
        fvm::laplacian(sigmaf, VRe)
//       - (gradAlpha&AIm)
    );

    VReEqn.setReference(refCell, refValue);

    // Imaginary part
    fvScalarMatrix VImEqn
    (
        fvm::laplacian(sigmaf, VIm)
//       + (gradAlpha&ARe)
    );

    VImEqn.setReference(refCell, refValue);

    // Remove non-conducting cells from solution of V
    // TODO: Use labelHashSet!
    label nVzeroCells = 0;
    forAll(mesh.cells(), cellI)
    {
        if
        (
            (mag(gradAlpha.internalField()[cellI]) < SMALL)
         && (sigma.internalField()[cellI] < SMALL)
        )
        {
            nVzeroCells++;
        }
    }
    labelList VzeroCells(nVzeroCells, -1);
    nVzeroCells = 0;
    forAll(mesh.cells(), cellI)
    {
        if
        (
            (mag(gradAlpha.internalField()[cellI]) < SMALL)
         && (sigma.internalField()[cellI] < SMALL)
        )
        {
            VzeroCells[nVzeroCells] = cellI;
            nVzeroCells++;
        }
    }
    scalarField VzeroCellValues(VzeroCells.size(), 0.0);
    VReEqn.setValues(VzeroCells, VzeroCellValues);
    VImEqn.setValues(VzeroCells, VzeroCellValues);
    
