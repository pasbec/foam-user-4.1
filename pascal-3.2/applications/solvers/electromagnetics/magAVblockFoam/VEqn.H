
//     // Reference to mesh data
//     const labelList& owner = mesh.owner();
//     const labelList& neighbour = mesh.neighbour();
//
//     // - Set an arbitrarily (high) value for sigmaf
//     //   for all faces which are in contact with
//     //   cells of zero (or very small) conductivity
//     //
//     // - Set alpha for all conductor faces to zero
//     //   in order to remove fluxes in VRe/VIm-equations
//     //   and achieve a zero-gradient boundary-condition
//     //   Set alpha for all faces between non-conducting
//     //   regions to zero
//     //
//     scalar sigmaNonConducting = max(mag(sigma)).value();
// //     scalar sigmaNonConducting = 1;
//     forAll(owner, faceI)
//     {
//         if
//         (
//             (sigma.internalField()[owner[faceI]] < VSMALL)
//          || (sigma.internalField()[neighbour[faceI]] < VSMALL)
//         )
//         {
//             sigmaf.internalField()[faceI] = sigmaNonConducting;
//             alphaf.internalField()[faceI] = 0;
//         }
//     }

    // Make sure sigmaf does not get zero
    scalar sigmafNonConducting = 1e-6;
    forAll(sigmaf.internalField(), faceI)
    {
        if
        (sigmaf.internalField()[faceI] < VSMALL)
        {
            sigmaf.internalField()[faceI] = sigmafNonConducting;
        }
    }



    // Real part

    fvScalarMatrix VReEqn
    (
      - fvm::laplacian(AVscale*sigmaf, VRe, "laplacian(sigmaf,VRe)")
//       + fvc::div(alpha*AIm)
    );

//     VReEqn.setReference(0, 0.0);

//     VReEqn.relax(0.9999);



    // Imaginary part
    fvScalarMatrix VImEqn
    (
      - fvm::laplacian(AVscale*sigmaf, VIm, "laplacian(sigmaf,VIm)")
//       - fvc::div(alpha*ARe)
    );

//     VImEqn.setReference(0, 0.0);

//     VImEqn.relax(0.9999);
