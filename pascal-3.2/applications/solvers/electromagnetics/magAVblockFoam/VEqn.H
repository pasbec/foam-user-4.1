
    // Calculate sigmaNoZero
    volScalarField sigmaNoZero
    (
        "sigmaNoZero",
        sigma
    );

    // Manipulate sigmaNoZero:
    //
    //  We will replace "zero"-values with very small
    //  values at all faces in order to avoid a singular
    //  solution
    {
        // Parameters
        scalar sigmaNoZeroRef = 100;

        // Manipulate internal faces
        forAll(sigmaNoZero.internalField(), cellI)
        {
            if
            (sigmaNoZero.internalField()[cellI] < sigmaNoZeroRef)
            {
                sigmaNoZero.internalField()[cellI] = sigmaNoZeroRef;
            }
        }

        // Manipulate boundary faces
        sigmaNoZero.boundaryField() == sigmaNoZeroRef;
    }

    // Interpolate sigma to face centers
    surfaceScalarField sigmaf
    (
        "sigmaf",
        fvc::interpolate(sigmaNoZero, "interpolate(sigma)")
    );

    // Calculate gradient of alpha
    volVectorField gradAlpha
    (
        "gradAlpha",
        omega * fvc::grad(sigma)
    );



    // Reference settings
    label refCell = mesh.findCell(point(0.00, 0.00, 0.01));
//     label refCell = mesh.findCell(point(0.034999, 0.00, 0.01));
//     label refCell = mesh.findCell(point(0.034999, 0.074999, 0.01));
//     label refCell = mesh.findCell(point(0.15, 0.11, 0.01));
    scalar refValue = 0.0;

//     Sout << "refCell = " << refCell << endl << endl;
//     Sout << "VRe.needReference() = " << VRe.needReference() << endl;
//     Sout << "VIm.needReference() = " << VIm.needReference() << endl;

    // Real part
    fvScalarMatrix VReEqn
    (
        fvm::laplacian(sigmaf, VRe, "laplacian(VRe)")
//       - fvc::div(alpha*AIm)
    );

    VReEqn.setReference(refCell, refValue);
//     VReEqn.relax(0.999);
// // TEST: Begin
//
//     if (VReEqn.psi().needReference())
//     {
//         Info << "VReEqn.source()[refCell] = " << VReEqn.source()[refCell] << endl;
//         Info << "VReEqn.diag()[refCell] = " << VReEqn.diag()[refCell] << endl;
//
// //         VReEqn.source()[refCell] += 1*VReEqn.diag()[refCell]*refValue;
// //         VReEqn.diag()[refCell] += 1*VReEqn.diag()[refCell];
//
//         VReEqn.source()[refCell] += 1.0*refValue;
//         VReEqn.diag()[refCell] += 1.0;
//
//         Info << "VReEqn.source()[refCell] = " << VReEqn.source()[refCell] << endl;
//         Info << "VReEqn.diag()[refCell] = " << VReEqn.diag()[refCell] << endl;
//     }
//
// // TEST: End

    // Imaginary part
    fvScalarMatrix VImEqn
    (
        fvm::laplacian(sigmaf, VIm, "laplacian(VIm)")
//       + fvc::div(alpha*ARe)
    );

    VImEqn.setReference(refCell, refValue);
//     VImEqn.relax(0.999);
// // TEST: Begin
//
//     if (VReEqn.psi().needReference())
//     {
//         Info << "VImEqn.source()[refCell] = " << VImEqn.source()[refCell] << endl;
//         Info << "VImEqn.diag()[refCell] = " << VImEqn.diag()[refCell] << endl;
//
// //         VImEqn.source()[refCell] += 1*VImEqn.diag()[refCell]*refValue;
// //         VImEqn.diag()[refCell] += 1*VImEqn.diag()[refCell];
//
//         VImEqn.source()[refCell] += 1.0*refValue;
//         VImEqn.diag()[refCell] += 1.0;
//
//         Info << "VImEqn.source()[refCell] = " << VImEqn.source()[refCell] << endl;
//         Info << "VImEqn.diag()[refCell] = " << VImEqn.diag()[refCell] << endl;
//     }
//
// // TEST: End

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

