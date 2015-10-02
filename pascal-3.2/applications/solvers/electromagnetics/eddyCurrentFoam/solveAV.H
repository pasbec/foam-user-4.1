
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
        // Interpolate A to face centers
        surfaceScalarField sARe
        (
            "sARe",
            fvc::interpolate(ARe) & mesh.Sf()
        );
        surfaceScalarField sAIm
        (
            "sAIm",
            fvc::interpolate(AIm) & mesh.Sf()
        );

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

            // Derived mesh geometry
            const vectorField Sn = Sf/magSf;

            forAll (cnPNfaces, facei)
            {
                label faceI = cnPNfaces[facei];

                // Calculate face centre/cell centre distance Pf
                scalar Pf = mag(Sn[faceI] & (Cf[faceI] - C[own[faceI]]));

                sVRe[faceI] = VRe[own[faceI]] + Pf/magSf[faceI] * omega.value() * sAIm[faceI];
                sVIm[faceI] = VIm[own[faceI]] - Pf/magSf[faceI] * omega.value() * sARe[faceI];
            }

            forAll (cnNPfaces, facei)
            {
                label faceI = cnNPfaces[facei];

                // Calculate face centre/cell centre distance fN
                scalar fN = mag(Sn[faceI] & (C[nei[faceI]] - Cf[faceI]));

                sVRe[faceI] = VRe[nei[faceI]] - fN/magSf[faceI] * omega.value() * sAIm[faceI];
                sVIm[faceI] = VIm[nei[faceI]] + fN/magSf[faceI] * omega.value() * sARe[faceI];
            }
        }

        // Update corrected gradient
        VReGrad = fvc::grad(sVRe);
        VImGrad = fvc::grad(sVIm);

        // Remove all gradient values from
        // non-conducting region
        VReGrad *= sigmaIc;
        VImGrad *= sigmaIc;
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
