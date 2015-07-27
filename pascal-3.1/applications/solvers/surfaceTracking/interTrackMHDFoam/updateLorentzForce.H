
    if (lorentzForceUpdateSwitch)
    {
        if (runTime.outputTime())
        {
            lorentzForceUpdateOutputTimeIndex += 1;
            lorentzForceTimeDict.set<int>("outputTimeIndex", lorentzForceUpdateOutputTimeIndex);

            if( lorentzForceUpdateOutputTimeIndex % lorentzForceUpdateOutputCycle == 0)
            {
#               include "updateLorentzForceRecalculate.H"
#               include "updateLorentzForcePressureCorr.H"

                lorentzForceTimeDict.set<scalar>("lastUpdateTime", runTime.value());
            }
        }
        else
        {
            lorentzForceUpdateTimeCycleLastUpdateTime = readScalar(lorentzForceTimeDict.lookup("lastUpdateTime"));

            if
            (
                (
                    runTime.timeIndex() % lorentzForceUpdateTimeStepCycle == 0
                )
                ||
                (
                    runTime.value() - lorentzForceUpdateTimeCycleLastUpdateTime > lorentzForceUpdateTimeCycle
                )
            )
            {
#               include "updateLorentzForceRecalculate.H"
#               include "updateLorentzForcePressureCorr.H"

                lorentzForceTimeDict.set<scalar>("lastUpdateTime", runTime.value());
            }
        }

        Info << nl;
    }