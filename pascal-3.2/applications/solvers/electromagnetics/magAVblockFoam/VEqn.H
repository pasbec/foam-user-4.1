
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
    //  First, we will replace zero values with very small
    //  values at all faces in order to avoid a singular
    //  solution
    //
    //  For all faces connecting cells with zero sigma
    //  (conductivity) and all boundary faces at infinity
    //  we now set the face interpolated data to some
    //  artificial (high) value. This will give an arbitrary
    //  solution for VRe/VIm in the region where sigma is
    //  actually zero and VRe/VIm. We need this in order to
    //  set some boundary conditions for VRe/VIm.
    {
        // Reference to mesh data
//         const unallocLabelList& own = mesh.owner();
//         const unallocLabelList& nei = mesh.neighbour();

        // Parameters
        scalar sigmafNonConducting = SMALL;
//         scalar sigmafDummy = max(sigmaf).value();

        // Manipulate internal faces
        forAll(sigmaf.internalField(), faceI)
        {
            if
            (sigmaf.internalField()[faceI] < SMALL)
            {
                sigmaf.internalField()[faceI] = sigmafNonConducting;

//                 label fOwn = own[faceI];
//                 label fNei = nei[faceI];
//
//                 if
//                 (
//                    (sigma.internalField()[fOwn] < SMALL)
//                 && (sigma.internalField()[fNei] < SMALL)
//                 && (mag(gradAlpha.internalField()[fOwn]) < SMALL)
//                 && (mag(gradAlpha.internalField()[fNei]) < SMALL)
//                 )
//                 {
//                     sigmaf.internalField()[faceI] = sigmafDummy;
//                 }
            }
        }

//         // Manipulate boundary faces
//         sigmaf.boundaryField() = sigmafDummy;
    }

    // Reference settings
    label refCell = mesh.findCell(point(0.00, 0.00, 0.01));
    scalar refValue = 0.0;

    // Real part
    fvScalarMatrix VReEqn
    (
        scaleV * fvm::laplacian(sigmaf, VRe)
//       - scaleV * (gradAlpha&AIm)
    );

    VReEqn.setReference(refCell, refValue);

    // Imaginary part
    fvScalarMatrix VImEqn
    (
        scaleV * fvm::laplacian(sigmaf, VIm)
//       + scaleV * (gradAlpha&ARe)
    );

    VImEqn.setReference(refCell, refValue);
