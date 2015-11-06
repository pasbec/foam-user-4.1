
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

            // Solve for A
#           include "solveAVsegregatedA.H"

            // Solve for V
#           include "solveAVsegregatedV.H"
        }
    }
    else if (AVmethod == "segregated")
    {
        int nCorr = AVmethodDict.lookupOrDefault<int>("corr",10);

        forAll (labelList(nCorr), corr)
        {
            // Assemble base-matrix for A
#           include "AEqn.H"

            // Solve for A
#           include "solveAVsegregatedA.H"

            // Assemble base-matrix for V
#           include "VEqn.H"

            // Solve for V
#           include "solveAVsegregatedVall.H"
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
#           include "VEqn.H"

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
