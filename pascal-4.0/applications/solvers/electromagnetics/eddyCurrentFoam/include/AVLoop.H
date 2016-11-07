
    {
        if (Control::debug)
        {
            Info<< Control::typeName << " | AVLoop.H : "
                << "Commencing EDDYCURRENT A-V loop."
                << endl;
        }

        // Init
#       include "AVLoop_sigmaUpdate.H"
#       include "AVLoop_AUpdate.H"
#       include "AVLoop_VGradUpdate.H"

        // AV loop
        while (control.loop())
        {
            while (control.subLoop())
            {
#               include "AVLoop_sigmaUpdate.H"

#               include "AVLoop_AEqn.H"

#               include "AVLoop_AUpdate.H"

#               include "AVLoop_VEqn.H"

#               include "AVLoop_VGradUpdate.H"
            }
        }
    }
