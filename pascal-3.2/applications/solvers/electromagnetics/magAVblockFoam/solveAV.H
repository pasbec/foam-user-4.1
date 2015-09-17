
    // Assemble matrix for A
#   include "AEqn.H"

    // Assemble matrix for V
#   include "VEqn.H"

    // Solve AV-system
    if (AVmethod == "subset")
    {
        int nSubsetCorr = AVmethodDict.lookupOrDefault<int>("corr",10);

        for
        (
            int subsetCorr=0;
            subsetCorr<nSubsetCorr;
            subsetCorr++
        )
        {
#           include "solveAVsubset.H"
        }
    }
    else if (AVmethod == "segregated")
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
    else if (AVmethod == "coupledSegregated")
    {
        int nCoupledCorr = AVmethodDict.lookupOrDefault<int>("coupledCorr",2);
        int nSegregatedCorr = AVmethodDict.lookupOrDefault<int>("segregatedCorr",5);

        for
        (
            int coupledCorr=0;
            coupledCorr<nCoupledCorr;
            coupledCorr++
        )
        {
#           include "solveAVcoupled.H"
        }

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
    else
    {
        FatalErrorIn("AVEqn.H")
        << "Method for solving AV-system '" << AVmethod
        << "' is not supported"
        << abort(FatalError);
    }
    
    jRe ==   alpha * AIm + alpha * fvc::grad(GIm) - sigma * fvc::grad(VRe);
    jIm == - alpha * ARe - alpha * fvc::grad(GRe) - sigma * fvc::grad(VIm);

    BRe == fvc::curl(ARe);
    BIm == fvc::curl(AIm);
