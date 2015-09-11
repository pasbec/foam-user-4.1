
    // Assemble matrix for A
#   include "AEqn.H"

    // Assemble matrix for V
#   include "VEqn.H"

    // Solve AV-system
    if (AVmethod == "coupled")
    {
        int nCoupledCorr = AVmethodDict.lookupOrDefault<int>("coupledCorr",4);

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
    else if (AVmethod == "segregated")
    {
        int nSegregatedCorr = AVmethodDict.lookupOrDefault<int>("segregatedCorr",10);

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
