
    // Get reference to indicator pointers
    volScalarField& indicatorVc = *indicatorVcPtr;
//     volScalarField& indicatorVn = *indicatorVnPtr;

    // AV iterations
    int nCorr = AVdict.lookupOrDefault<int>("corr",4);

    forAll (labelList(nCorr), corr)
    {
        // Solve for A
#       include "AEqn.H"

        // Solve for V
#       include "VEqn.H"
    }
