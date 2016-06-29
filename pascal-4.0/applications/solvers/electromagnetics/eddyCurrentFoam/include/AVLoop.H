
    // Init
#   include "AVLoop_AUpdate.H"
#   include "AVLoop_VGradUpdate.H"

    // AV loop
    while (control.loop())
    {
        while (control.subLoop())
        {
#           include "AVLoop_AEqn.H"

#           include "AVLoop_AUpdate.H"

// TODO: Fix relaxation
//             if (control.meshIs3D())
//             {
//                 // Relax A
//                 control.relax(globalStorage.ARe());
//                 control.relax(globalStorage.AIm());
//             }

#           include "AVLoop_VEqn.H"

#           include "AVLoop_VGradUpdate.H"
        }
    }
