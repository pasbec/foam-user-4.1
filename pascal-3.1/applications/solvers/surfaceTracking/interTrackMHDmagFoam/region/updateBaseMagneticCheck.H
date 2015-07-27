
// TODO [High]: Magnetic update does not work correctly and happens
//              way to often, especially before first write

    if
    (
        mfUpdateSwitch
     || mfUpdateCounter == 0
    )
    {
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
                mfUpdateCounter == 0
            )
            ||
            (
                (mfUpdateOutputTimeIndex
                    % mfUpdateOutputCycle) == 0
            )
            ||
            (
                (runTime.timeIndex()
                    % mfUpdateTimeStepCycle) == 0
            )
            ||
            (
                (runTime.value()-mfUpdateLastUpdateTime)
                    - mfUpdateTimeCycle > 0
            )
        )
        {
            mfUpdateCounter += 1;
            mfTimeDict.set<int>
            (
                "updateCounter",
                mfUpdateCounter
            );

            mfUpdateLastUpdateTime = runTime.value();
            mfTimeDict.set<scalar>
            (
                "lastUpdateTime",
                mfUpdateLastUpdateTime
            );

            mfUpdate = true;
        }
        else
        {
            mfUpdate = false;
        }
    }
