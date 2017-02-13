
    // Init
#   include "AVLoop_sigmaUpdate.H"
#   include "AVLoop_murUpdate.H"

// TEST
    {
        SM_REGIONSCOPE(DEFAULT);

        for (int AEqnI = 0; AEqnI < 2; AEqnI++)
        {
            volVectorField* APtr = NULL;
            surfaceVectorField* gradNPhiPtr = NULL;

            if (AEqnI == 0)
            {
                APtr = &globalStorage.ARe()[DEFAULT];
                gradNPhiPtr = &storage.gradAReNPhi();
            }
            else
            {
                APtr = &globalStorage.AIm()[DEFAULT];
                gradNPhiPtr = &storage.gradAImNPhi();
            }

            volVectorField& A = *APtr;
            surfaceVectorField& gradNPhi = *gradNPhiPtr;

            surfaceTensorField gradNf =
                fvc::interpolate
                (
                    fvc::grad(linearInterpolate(A)),
                    "interpolate(grad(" + A.name() + "))"
                );

            gradNPhi = -gradNf & mesh.Sf();

            gradNPhi += tr(gradNf) * mesh.Sf();
        }
    }

#   include "AVLoop_AUpdate.H"
#   include "AVLoop_VGradInit.H"
#   include "AVLoop_VGradUpdate.H"
