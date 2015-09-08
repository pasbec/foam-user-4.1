
    // Make sure sigmaf does not get zero
    scalar sigmafNonConducting = SMALL;
    forAll(sigmaf.internalField(), faceI)
    {
        if
        (sigmaf.internalField()[faceI] < SMALL)
        {
            sigmaf.internalField()[faceI] = sigmafNonConducting;
        }
    }

    // Reference settings
    label refCell = mesh.findCell(point(0.00, 0.00, 0.01));
    scalar refValue = 0.0;

    // Real part

    fvScalarMatrix VReEqn
    (
        Vscale * fvm::laplacian(sigmaf, VRe, "laplacian(sigmaf,VRe)")
//       - Vscale * (gradAlpha&AIm)
    );

//     VReEqn.setReference(refCell, refValue);



    // Imaginary part
    fvScalarMatrix VImEqn
    (
        Vscale * fvm::laplacian(sigmaf, VIm, "laplacian(sigmaf,VIm)")
//       + Vscale * (gradAlpha&ARe)
    );

//     VImEqn.setReference(refCell, refValue);
