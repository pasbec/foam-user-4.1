
// TODO [High]: Magnetic update does not work correctly and happens
//              way to often, especially before first write

    Info << "[DEBUG] runTime.value() = " << runTime.value() << endl;
    Info << "[DEBUG] runTime.timeIndex() = " << runTime.timeIndex() << endl;
    Info << "[DEBUG] runTime.outputTime() = " << runTime.outputTime() << endl;

    Info << "[DEBUG] mfUpdateSwitch = " << mfUpdateSwitch << endl;
    Info << "[DEBUG] mfUpdateOutputCycle = " << mfUpdateOutputCycle << endl;
    Info << "[DEBUG] mfUpdateTimeStepCycle = " << mfUpdateTimeStepCycle << endl;
    Info << "[DEBUG] mfUpdateTimeCycle = " << mfUpdateTimeCycle << endl;

    Info << "[DEBUG] mfUpdateCounter = " << mfUpdateCounter << endl;
    Info << "[DEBUG] mfUpdateLastUpdateTime = " << mfUpdateLastUpdateTime << endl;
    Info << "[DEBUG] mfUpdateOutputTimeIndex = " << mfUpdateOutputTimeIndex << endl;

    if (mfUpdateSwitch)
    {
        // Increase output index if due and write to dictionary
        if (runTime.outputTime())
        {
            mfUpdateOutputTimeIndex += 1;
            mfTimeDict.set<int>
            (
                "outputTimeIndex",
                mfUpdateOutputTimeIndex
            );
        }

        Info << "[DEBUG] mfUpdateOutputTimeIndex = " << mfUpdateOutputTimeIndex << endl;
        Info << "[DEBUG] mfUpdateCounter = " << mfUpdateCounter << endl;
        Info << "[DEBUG] (mfUpdateOutputTimeIndex % mfUpdateOutputCycle) = " << (mfUpdateOutputTimeIndex % mfUpdateOutputCycle) << endl;
        Info << "[DEBUG] (runTime.timeIndex() % mfUpdateTimeStepCycle) = " << (runTime.timeIndex() % mfUpdateTimeStepCycle) << endl;
        Info << "[DEBUG] (runTime.value() - mfUpdateLastUpdateTime - mfUpdateTimeCycle = " << (runTime.value()-mfUpdateLastUpdateTime - mfUpdateTimeCycle) << endl;

        if
        (
            (
                // Update has never run
                mfUpdateCounter == 0
            )
            ||
            (
                // Update triggered by output cycle interval is due
                (
                    mfUpdateOutputTimeIndex > 0
                )
                &&
                (
                    mfUpdateOutputTimeIndex % mfUpdateOutputCycle == 0
                )
            )
            ||
            (
                // Update triggered by time step cycle interval is due
                (
                    runTime.timeIndex() > 0
                )
                &&
                (
                    runTime.timeIndex() % mfUpdateTimeStepCycle == 0
                )
            )
            ||
            (
                // Update triggered by time cycle interval is due
                runTime.value() - mfUpdateLastUpdateTime > mfUpdateTimeCycle
            )
        )
        {
            mfUpdate = true;

            // Increase update counter and write to dictionary
            mfUpdateCounter += 1;
            mfTimeDict.set<int>
            (
                "updateCounter",
                mfUpdateCounter
            );

            // Update last update time and write to dictionary
            mfUpdateLastUpdateTime = runTime.value();
            mfTimeDict.set<scalar>
            (
                "lastUpdateTime",
                mfUpdateLastUpdateTime
            );
        }
        else
        {
            mfUpdate = false;
        }
    }

    Info << "[DEBUG] mfUpdate = " << mfUpdate << endl;
