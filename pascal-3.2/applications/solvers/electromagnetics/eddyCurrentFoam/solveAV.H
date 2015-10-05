
    // Get reference to indicator pointers
    volScalarField& indicatorVc = *indicatorVcPtr;
//     volScalarField& indicatorVn = *indicatorVnPtr;

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

        // Correct gradient only at the side of
        // conducting region
        {
            // Mesh cell-face addressing
            const unallocLabelList& own = mesh.owner();
            const unallocLabelList& nei = mesh.neighbour();

            // Mesh geometry
            const vectorField& Sf = mesh.Sf();
            const scalarField& magSf = mesh.magSf();
            const scalarField& w = mesh.weights();
            const scalarField& d = mesh.deltaCoeffs();

            // Derived mesh geometry
            const vectorField Sn = Sf/magSf;

            forAll (cnPNfaces, facei)
            {
                label faceI = cnPNfaces[facei];

                // Calculate face centre to cell centre distance Pf
                scalar Pf = (1.0 - w[faceI]) / d[faceI];

                // Calculate weight coeffs
                scalar wP = omega.value() * w[faceI];
                scalar wN = omega.value() - wP;

                sVRe[faceI] = VRe[own[faceI]];
                sVRe[faceI] += Pf * wP * (Sn[faceI] & AIm[own[faceI]]);
                sVRe[faceI] += Pf * wN * (Sn[faceI] & AIm[nei[faceI]]);

                sVIm[faceI] = VIm[own[faceI]];
                sVIm[faceI] -= Pf * wP * (Sn[faceI] & ARe[own[faceI]]);
                sVIm[faceI] -= Pf * wN * (Sn[faceI] & ARe[nei[faceI]]);
            }

            forAll (cnNPfaces, facei)
            {
                label faceI = cnNPfaces[facei];

                // Calculate face centre to cell centre distance fN
                scalar fN = w[faceI] / d[faceI];

                // Calculate weight coeffs
                scalar wP = omega.value() * w[faceI];
                scalar wN = omega.value() - wP;

                sVRe[faceI] = VRe[nei[faceI]];
                sVRe[faceI] -= fN * wP * (Sn[faceI] & AIm[own[faceI]]);
                sVRe[faceI] -= fN * wN * (Sn[faceI] & AIm[nei[faceI]]);

                sVIm[faceI] = VIm[nei[faceI]];
                sVIm[faceI] += fN * wP * (Sn[faceI] & ARe[own[faceI]]);
                sVIm[faceI] += fN * wN * (Sn[faceI] & ARe[nei[faceI]]);
            }
        }

        // Store old gradient
        VReGrad.storePrevIter();
        VImGrad.storePrevIter();

        // Update corrected gradient
        VReGrad = fvc::grad(sVRe) * indicatorVc;
        VImGrad = fvc::grad(sVIm) * indicatorVc;
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
