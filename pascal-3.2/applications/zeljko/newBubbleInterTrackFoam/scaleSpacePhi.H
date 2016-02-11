if(runTime.timeIndex() > 1)
{
    fvsPatchScalarField& spacePhi = phi.boundaryField()[spacePatchID];

    scalar inletFlux = gSum(neg(spacePhi)*spacePhi);

//     Pout << "inletFlux: " << inletFlux << endl;

    scalar outletFlux = gSum(pos(spacePhi)*spacePhi);

//     Pout << "outletFlux: " << outletFlux << endl;

    if(outletFlux < VSMALL)
    {
        outletFlux = VSMALL;
    }

    scalar outflowScaling = -inletFlux/outletFlux;

//     Pout << "outflowScaling: " << outflowScaling << endl;

    spacePhi += pos(spacePhi)*spacePhi*(outflowScaling - 1.0);

    Pout << "Space patch continuity: " 
        << gSum(phi.boundaryField()[spacePatchID]) << endl;

//     U.boundaryField()[spacePatchID] ==
//         U.boundaryField()[spacePatchID]
//       + pos(spacePhi)*U.boundaryField()[spacePatchID]*(outflowScaling - 1.0);
}
