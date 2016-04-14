
    // Init gradient of V in first iteration
    // in case of a three-dimensional mesh
    if (control.meshIs3D())
    {
#       include "AVLoop_VGradUpdate.H"
    }

    // AV loop
    while (control.loop())
    {
        while (control.subLoop())
        {
#           include "AVLoop_AEqn.H"

            // The electric scalar potential is only necessary
            // in case of a three-dimensional mesh
            if (control.meshIs3D())
            {
#               include "AVLoop_AUpdate.H"

// TODO: Fix relaxation
//                 // Relax A
//                 control.relax(globalStorage.ARe());
//                 control.relax(globalStorage.AIm());

#               include "AVLoop_VEqn.H"

#               include "AVLoop_VGradUpdate.H"
            }
        }
    }
