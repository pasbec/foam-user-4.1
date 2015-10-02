
    // Assemble matrix for A
#   include "AEqn.H"

    // Assemble matrix for V
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

    // Update corrected gradient of V
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

                // Linear interpolate A
                vector sAReSn = wP * ARe[own[faceI]] + wN * ARe[nei[faceI]];
                vector sAImSn = wP * AIm[own[faceI]] + wN * AIm[nei[faceI]];

                sVRe[faceI] = VRe[own[faceI]] + omega.value() * Pf * (Sn[faceI] & sAImSn);
                sVIm[faceI] = VIm[own[faceI]] - omega.value() * Pf * (Sn[faceI] & sAReSn);
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

                // Linear interpolate A
                vector sAReSn = wP * ARe[own[faceI]] + wN * ARe[nei[faceI]];
                vector sAImSn = wP * AIm[own[faceI]] + wN * AIm[nei[faceI]];

                sVRe[faceI] = VRe[nei[faceI]] - omega.value() * fN * (Sn[faceI] & sAImSn);
                sVIm[faceI] = VIm[nei[faceI]] + omega.value() * fN * (Sn[faceI] & sAReSn);
            }
        }

        // Update corrected gradient
        VReGrad = fvc::grad(sVRe);
        VImGrad = fvc::grad(sVIm);

        // Remove all gradient values from
        // non-conducting region
        VReGrad *= indicatorVc;
        VImGrad *= indicatorVc;
    }



    // Calculate derived fields
    {
        // Magnetic field density
        BRe == fvc::curl(ARe);
        BIm == fvc::curl(AIm);

        // Eddy current density
        jRe ==   alpha * AIm - sigma * VReGrad;
        jIm == - alpha * ARe - sigma * VImGrad;
    }
