
    // Get reference to indicator pointers
    volScalarField& indicatorVc = *indicatorVcPtr;
    volScalarField& indicatorVn = *indicatorVnPtr;

    // Update corrected gradient of V
    // TODO: Put this somewhere else
    {
        // Interpolate V to face centers
        surfaceScalarField sVRe
        (
            "sVRe",
            fvc::interpolate(VRe)
        );
        surfaceScalarField sVIm
        (
            "sVIm",
            fvc::interpolate(VIm)
        );

        // Correct gradient
        {
            // Mesh cell-face addressing
            const unallocLabelList& own = mesh.owner();
            const unallocLabelList& nei = mesh.neighbour();

            // Mesh geometry
            const vectorField& Cf = mesh.Cf();
            const vectorField& C = mesh.C();
            const vectorField& Sf = mesh.Sf();
            const scalarField& magSf = mesh.magSf();
            const vectorField Sn = Sf/magSf;

            forAll (cnPNfaces, facei)
            {
                label faceI = cnPNfaces[facei];

                // Calculate face centre to cell centre distances Pf and fN
                scalar Pf = mag(Sn[faceI] & (Cf[faceI] - C[own[faceI]]));
                scalar fN = mag(Sn[faceI] & (C[nei[faceI]] - Cf[faceI]));

                // Calculate weights
                scalar wP = fN / (Pf+fN);
                scalar wN = Pf / (Pf+fN);

                // Calculate gradient coeffs
                scalar cP = omega.value() * Pf * wP;
                scalar cN = omega.value() * Pf * wN;

                sVRe[faceI] = VRe[own[faceI]] + cP * (Sn[faceI] & AIm[own[faceI]])
                                              + cN * (Sn[faceI] & AIm[nei[faceI]]);
                sVIm[faceI] = VIm[own[faceI]] - cP * (Sn[faceI] & ARe[own[faceI]])
                                              - cN * (Sn[faceI] & ARe[nei[faceI]]);
            }

            forAll (cnNPfaces, facei)
            {
                label faceI = cnNPfaces[facei];

                // Calculate face centre to cell centre distances Pf and fN
                scalar Pf = mag(Sn[faceI] & (Cf[faceI] - C[own[faceI]]));
                scalar fN = mag(Sn[faceI] & (C[nei[faceI]] - Cf[faceI]));

                // Calculate weights
                scalar wP = fN / (Pf+fN);
                scalar wN = Pf / (Pf+fN);

                // Calculate gradient coeffs
                scalar cP = omega.value() * fN * wP;
                scalar cN = omega.value() * fN * wN;

                sVRe[faceI] = VRe[nei[faceI]] - cP * (Sn[faceI] & AIm[own[faceI]])
                                              - cN * (Sn[faceI] & AIm[nei[faceI]]);
                sVIm[faceI] = VIm[nei[faceI]] + cP * (Sn[faceI] & ARe[own[faceI]])
                                              + cN * (Sn[faceI] & ARe[nei[faceI]]);
            }
        }

        // Store old gradient
        VReGrad.storePrevIter();
        VImGrad.storePrevIter();

        // Update corrected gradient
        VReGrad = fvc::grad(sVRe);
        VImGrad = fvc::grad(sVIm);

        // Remove all gradient values from
        // non-conducting region
        VReGrad *= indicatorVc;
        VImGrad *= indicatorVc;
    }

    // Assemble base-matrix for A
#   include "AEqn.H"

    // Assemble base-matrix for V
#   include "VEqn.H"

    // Solve AV-system
    if (AVmethod == "segregated")
    {
        int nSegregatedCorr = AVmethodDict.lookupOrDefault<int>("corr",10);

        for
        (
            int segregatedCorr=0;
            segregatedCorr<nSegregatedCorr;
            segregatedCorr++
        )
        {
#           include "solveAVsegregated.H"
        }
    }
    else if (AVmethod == "coupled")
    {
        int nCoupledCorr = AVmethodDict.lookupOrDefault<int>("corr",4);

        for
        (
            int coupledCorr=0;
            coupledCorr<nCoupledCorr;
            coupledCorr++
        )
        {
#           include "solveAVcoupled.H"
        }
    }
    else
    {
        FatalErrorIn("AVEqn.H")
        << "Method for solving AV-system '" << AVmethod
        << "' is not supported"
        << abort(FatalError);
    }

    // Calculate derived fields
    {
        // Magnetic field density
        BRe == fvc::curl(ARe);
        BIm == fvc::curl(AIm);

        // Eddy current density
        jRe ==   alpha * AIm - sigma * VReGrad;
        jIm == - alpha * ARe - sigma * VImGrad;

        // TODO: Remove after debug
        AReDiv == fvc::div(ARe);
        AImDiv == fvc::div(AIm);
    }
