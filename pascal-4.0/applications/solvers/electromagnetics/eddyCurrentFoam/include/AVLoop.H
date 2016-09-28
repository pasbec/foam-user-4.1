
    {
        if (Control::debug)
        {
            Info << "eddyCurrentApp::Control | UpLoop.H : "
                << "Commencing EDDYCURRENT A-V loop."
                << endl;
        }

        // Init
#       include "AVLoop_AUpdate.H"
#       include "AVLoop_VGradUpdate.H"

        // AV loop
        while (control.loop())
        {
            while (control.subLoop())
            {
#               include "AVLoop_AEqn.H"

#               include "AVLoop_AUpdate.H"

#               include "AVLoop_VEqn.H"

#               include "AVLoop_VGradUpdate.H"
            }
        }
    }
