
#   include "volContinuity.H"

    // Surface continuity error
// TODO [Temp]
//     scalarField netPhiA = interface.continuityErrorA();
    scalarField netPhiA =
        phi.boundaryField()[interface.aPatchID()]
      - fvc::meshPhi(U)().boundaryField()[interface.aPatchID()];

    Info << "Surface continuity error: " << gSum(netPhiA) << endl;

    // Relative volume error
    scalar relVolErr = mag( gSum(netPhiA) / gSum(mesh.cellVolumes()) );

    // Print relative volume error
    Info << "Relative volume error: " << relVolErr << endl;

    // Prepare PIMPLE loop break if conditions met
    if
    (
        (
            !outerCorrBreak
        )
        &&
        (
            outerCorr + 2 > nMinOuterCorr
        )
        &&
        (
            relVolErr <= maxRelVolErr
        )
    )
    {
        // Break info
        Info << "Target for relative volume error reached." << endl
             << "Exiting outer PIMPLE loop after next iteration." << endl;

        // Prepare break
        nOuterCorr = outerCorr + 2;
        outerCorrBreak = true;
    }

    Info << endl;