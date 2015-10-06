
    // Get reference to indicator pointers
    volScalarField& indicatorVc = *indicatorVcPtr;
//     volScalarField& indicatorVn = *indicatorVnPtr;

    // Update gradient of V
    {
        // Store old gradient
        VReGrad.storePrevIter();
        VImGrad.storePrevIter();

        VReGrad = fvc::grad(VRe);
        VImGrad = fvc::grad(VIm);
    }

    // Assemble base-matrix for A
#   include "AEqn.H"

    // Assemble base-matrix for V
#   include "VEqn.H"

    // Solve AV-system
    if (AVmethod == "segregated")
    {
        int nCorr = AVmethodDict.lookupOrDefault<int>("corr",10);

        forAll (labelList(nCorr), corr)
        {
#           include "solveAVsegregated.H"
        }
    }
    else if (AVmethod == "coupled")
    {
        int nCorr = AVmethodDict.lookupOrDefault<int>("corr",4);

        forAll (labelList(nCorr), corr)
        {
#           include "solveAVcoupled.H"
        }
    }
    else
    {
        FatalErrorIn("AVEqn.H")
        << "Method for solving AV-system '"
        << AVmethod
        << "' is not supported!"
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
