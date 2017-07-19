{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phic = mag(phi/mesh.magSf());
    phic = min(interface.cAlpha()*phic, max(phic));
    surfaceScalarField phir = phic*interface.nHatf();

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        surfaceScalarField phiAlpha =
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
                alpha1,
                alpharScheme
            );

        MULES::explicitSolve(alpha1, phi, phiAlpha, 1, 0);

// // // TODO: Test
        dictionary& pisoDict = mesh.solutionDict().subDict("PISO");
        if(pisoDict.lookupOrDefault("cutCompression", false))
        {
            volScalarField magGradAlpha("magGradAlpha", mag(fvc::grad(alpha1)));
            magGradAlpha /= max(magGradAlpha)
                + dimensionedScalar(word(), magGradAlpha.dimensions(), VSMALL);

            scalar cAlphaCut = pisoDict.lookupOrDefault("cAlphaCut", 0.0);
            scalar cAlphaCutRelErr = pisoDict.lookupOrDefault("cAlphaCutRelErr", 1e-4);

Info<< "DEBUG | cAlpha = " << interface.cAlpha() << endl;
Info<< "DEBUG | cAlphaCut = " << cAlphaCut << endl;
Info<< "DEBUG | cAlphaCutRelErr = " << cAlphaCutRelErr << endl;

            scalar scaleCut = cAlphaCut;

            volScalarField scaleCutAlpha1("scaleCutAlpha1", (1.0 + scaleCut*(1.0-magGradAlpha))*2.0*(alpha1 - 0.5));
            scaleCutAlpha1 = 0.5*min(max(scaleCutAlpha1, scalar(-1)), scalar(1)) + 0.5;
            scaleCutAlpha1.correctBoundaryConditions();
            dimensionedScalar scaleCutAlpha1Vfrac = scaleCutAlpha1.weightedAverage(mesh.V());
            dimensionedScalar relErrScaleCutAlpha1Vfrac = (scaleCutAlpha1Vfrac-alpha1Vfrac)/alpha1Vfrac;
// scaleCutAlpha1.write();

Info<< "DEBUG | scaleCutAlpha1Vfrac.value() = " << scaleCutAlpha1Vfrac.value() << endl;
Info<< "DEBUG | relErrScaleCutAlpha1Vfrac.value() = " << relErrScaleCutAlpha1Vfrac.value() << endl;

            volScalarField shiftCutAlpha1("shiftCutAlpha1", scaleCutAlpha1);

            if (mag(relErrScaleCutAlpha1Vfrac).value() > cAlphaCutRelErr)
            {
                scalar  shiftCut = cAlphaCut;

                shiftCutAlpha1 *= (1.0 + shiftCut*(1.0-magGradAlpha));

                if (relErrScaleCutAlpha1Vfrac.value() > 0.0)
                {
                    shiftCutAlpha1 -= shiftCut*(1.0-magGradAlpha);
                }

                shiftCutAlpha1 = min(max(shiftCutAlpha1, scalar(0)), scalar(1));
                shiftCutAlpha1.correctBoundaryConditions();
                dimensionedScalar shiftCutAlpha1Vfrac = shiftCutAlpha1.weightedAverage(mesh.V());
                dimensionedScalar relErrShiftCutAlpha1Vfrac = (shiftCutAlpha1Vfrac-alpha1Vfrac)/alpha1Vfrac;
// shiftCutAlpha1.write();

Info<< "DEBUG | shiftCutAlpha1Vfrac.value() = " << shiftCutAlpha1Vfrac.value() << endl;
Info<< "DEBUG | relErrShiftCutAlpha1Vfrac.value() = " << relErrShiftCutAlpha1Vfrac.value() << endl;
            }

            alpha1 = shiftCutAlpha1;
            alpha1.correctBoundaryConditions();
        }

        rhoPhi = phiAlpha*(transport.rho1() - transport.rho2()) + phi*transport.rho2();
    }

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;
}
