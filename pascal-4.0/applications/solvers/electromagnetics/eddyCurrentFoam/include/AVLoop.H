
    {
        if (Control::debug)
        {
            Info<< Control::typeName << " | AVLoop.H : "
                << "Commencing EDDYCURRENT A-V loop."
                << endl;
        }

// TODO: Residual controlled loop size
        for (int corrI = 0; corrI < 20; corrI++)
        {
#           include "AVLoop_A0Eqn.H"

#           include "AVLoop_A0Update.H"
        }

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
