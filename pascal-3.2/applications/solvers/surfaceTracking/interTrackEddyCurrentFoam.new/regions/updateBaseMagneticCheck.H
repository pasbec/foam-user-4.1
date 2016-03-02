
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
