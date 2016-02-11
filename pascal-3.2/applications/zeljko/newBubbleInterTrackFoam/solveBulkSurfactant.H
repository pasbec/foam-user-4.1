
if(CPtr)
{
    volScalarField& C = *CPtr;

    const dimensionedScalar& D 
        = interface.surfactant().surfactBulkDiffusion();

    scalar ka = interface.surfactant().surfactAdsorptionCoeff().value();
    scalar kb = interface.surfactant().surfactDesorptionCoeff().value();
    scalar CsInf = interface.surfactant().surfactSaturatedConc().value();

    scalarField& Cs = 
        interface.surfactantConcentration().internalField();

    scalarField CP = 
        C.boundaryField()[interface.aPatchID()].patchInternalField();

    const scalarField& Cfs = C.boundaryField()[interface.aPatchID()];

    scalarField dn = 1.0/mesh.boundary()[interface.aPatchID()].deltaCoeffs();

    if
    (
        C.boundaryField()[interface.aPatchID()].type() 
     == fixedGradientFvPatchScalarField::typeName
    )
    {
        fixedGradientFvPatchScalarField& CA =
            refCast<fixedGradientFvPatchScalarField>
            (
                C.boundaryField()[interface.aPatchID()]
            );

        CA.gradient() = (ka*kb*Cs - ka*(CsInf-Cs)*Cfs)/D.value();

//         CA.gradient() = 
//             (ka*kb*Cs - ka*(CsInf-Cs)*CP)
//            /(D.value() + dn*ka*(CsInf-Cs));
    }
//     else if
//     (
//         C.boundaryField()[interface.aPatchID()].type() 
//      == fixedGradientCorrectedFvPatchScalarField::typeName
//     )
//     {
//         fixedGradientCorrectedFvPatchScalarField& CA =
//             refCast<fixedGradientCorrectedFvPatchScalarField>
//             (
//                 C.boundaryField()[interface.aPatchID()]
//             );

//         CA.gradient() = (ka*kb*Cs - ka*(CsInf-Cs)*Cfs)/D.value();

// //         CA.gradient() = 
// //             (ka*kb*Cs - ka*(CsInf-Cs)*CP)
// //            /(D.value() + dn*ka*(CsInf-Cs));
//     }
    else
    {
        FatalErrorIn(args.executable())
            << "Bulk concentration boundary condition "
                << "at the freeSurface patch is not "
                << fixedGradientFvPatchScalarField::typeName
//                 << " or "
//                 << fixedGradientCorrectedFvPatchScalarField::typeName
                << exit(FatalError);
    }
    
    solve
    (
        fvm::ddt(C)
      + fvm::div(phi - fvc::meshPhi(rho,U), C, "div(phi,C)")
      - fvm::laplacian(D, C, "laplacian(D,C)")
    );

    CP = C.boundaryField()[interface.aPatchID()].patchInternalField();

    Info << gMax(CP) << ", "<< gAverage(CP) << endl;
}
