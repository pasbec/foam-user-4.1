// if(Pstream::master())
// {
    Info<< "Surface flux: sum local = "
        << gSum(mag(phi.boundaryField()[interface.aPatchID()]))
        << ", global = "
        << gSum(phi.boundaryField()[interface.aPatchID()]) << endl;

    scalarField netPhiA =
        phi.boundaryField()[interface.aPatchID()]
      - fvc::meshPhi(rho,U)().boundaryField()[interface.aPatchID()];

    Info<< "Surface continuity error : sum local = "
        << gSum(mag(netPhiA)) << ", global = " << gSum(netPhiA)
        << endl << endl;
// }
