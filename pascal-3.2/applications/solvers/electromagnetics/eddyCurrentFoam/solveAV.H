
    // Get reference to indicator pointers
    volScalarField& indicatorVc = *indicatorVcPtr;
//     volScalarField& indicatorVn = *indicatorVnPtr;

    // Solve AV-system
    if (AVmethod == "subset")
    {
        int nCorr = AVmethodDict.lookupOrDefault<int>("corr",4);

        forAll (labelList(nCorr), corr)
        {
            // Assemble base-matrix for A
#           include "AEqn.H"

#           include "solveAVsubset.H"
        }
    }
    else if (AVmethod == "segregated")
    {
        int nCorr = AVmethodDict.lookupOrDefault<int>("corr",10);

        forAll (labelList(nCorr), corr)
        {
            // Assemble base-matrix for A
#           include "AEqn.H"

            // Assemble base-matrix for V
#           include "VEqn.H"

#           include "solveAVsegregated.H"
        }
    }
    else if (AVmethod == "coupled")
    {
        int nCorr = AVmethodDict.lookupOrDefault<int>("corr",4);

        forAll (labelList(nCorr), corr)
        {
            // Assemble base-matrix for A
#           include "AEqn.H"

            // Assemble base-matrix for V
        #   include "VEqn.H"

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
