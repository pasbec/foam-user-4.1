
    // Calculate gradient of alpha
    volVectorField gradAlpha
    (
        "gradAlpha",
        omega * fvc::grad(sigma)
    );

    // Interpolate sigma to face centers
    surfaceScalarField sigmaf
    (
        "sigmaf",
        fvc::interpolate(sigma)
    );

    // Manipulate face interpolated sigma:
    //
    //  We will replace "zero"-values with very small
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
//     label refCell = mesh.findCell(point(0.00, 0.00, 0.01));
//     label refCell = mesh.findCell(point(0.034999, 0.00, 0.01));
    label refCell = mesh.findCell(point(0.034999, 0.074999, 0.01));
//     label refCell = mesh.findCell(point(0.15, 0.11, 0.01));
    scalar refValue = 0.0;

    // Real part
    fvScalarMatrix VReEqn
    (
        fvm::laplacian(sigmaf, VRe)
//       - (gradAlpha&AIm)
    );

//     VReEqn.setReference(refCell, refValue);
// TEST: Begin
    if (VReEqn.psi().needReference())
    {
        Info << "VReEqn.source()[refCell] = " << VReEqn.source()[refCell] << endl;
        Info << "VReEqn.diag()[refCell] = " << VReEqn.diag()[refCell] << endl;

        VReEqn.source()[refCell] += 1*VReEqn.diag()[refCell]*refValue;
        VReEqn.diag()[refCell] += 1*VReEqn.diag()[refCell];

//         VReEqn.source()[refCell] += 1.0*refValue;
//         VReEqn.diag()[refCell] += 1.0;

        Info << "VReEqn.source()[refCell] = " << VReEqn.source()[refCell] << endl;
        Info << "VReEqn.diag()[refCell] = " << VReEqn.diag()[refCell] << endl;
    }
// TEST: End

    // Imaginary part
    fvScalarMatrix VImEqn
    (
        fvm::laplacian(sigmaf, VIm)
//       + (gradAlpha&ARe)
    );

//     VImEqn.setReference(refCell, refValue);
// TEST: Begin
    if (VReEqn.psi().needReference())
    {
        Info << "VImEqn.source()[refCell] = " << VImEqn.source()[refCell] << endl;
        Info << "VImEqn.diag()[refCell] = " << VImEqn.diag()[refCell] << endl;

        VImEqn.source()[refCell] += 1*VImEqn.diag()[refCell]*refValue;
        VImEqn.diag()[refCell] += 1*VImEqn.diag()[refCell];

//         VImEqn.source()[refCell] += 1.0*refValue;
//         VImEqn.diag()[refCell] += 1.0;

        Info << "VImEqn.source()[refCell] = " << VImEqn.source()[refCell] << endl;
        Info << "VImEqn.diag()[refCell] = " << VImEqn.diag()[refCell] << endl;
    }
// TEST: End

// // TEST: Begin
//
//     // Create set for non-conducting cells
//     labelHashSet nonConductingCellsHashSet;
//     forAll(mesh.cells(), cellI)
//     {
//         if
//         (
//             (sigma.internalField()[cellI] < SMALL)
//          && (mag(gradAlpha.internalField()[cellI]) < SMALL)
//         )
//         {
//             nonConductingCellsHashSet.insert(cellI);
//         }
//     }
//     labelList nonConductingCells = nonConductingCellsHashSet.toc();
//
//     // Set values for non-conducting cells
//     scalarField nonConductingCellValues(nonConductingCells.size(), 0.0);
//     VReEqn.setValues(nonConductingCells, nonConductingCellValues);
//     VImEqn.setValues(nonConductingCells, nonConductingCellValues);
//
// // TEST: End

