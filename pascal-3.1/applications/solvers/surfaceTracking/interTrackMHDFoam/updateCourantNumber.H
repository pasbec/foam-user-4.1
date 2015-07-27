
#   include "CourantNo.H"
#   include "meshCourantNo.H"

    scalar convectiveCoNum = max(CoNum, meshCoNum);
    scalar interfaceCoNum = interface.maxCourantNumber();

    Info << "Courant Number of tracked surface max: " << interfaceCoNum << endl;

    scalar restrictConvectiveCoNum = convectiveCoNum/maxConvectiveCo;
    scalar restrictInterfaceCoNum = interfaceCoNum/maxInterfaceCo;

    if (restrictConvectiveCoNum < restrictInterfaceCoNum)
    {
        CoNum = interfaceCoNum;
        maxCo = min(maxInterfaceCo,interfaceCoNum/(restrictConvectiveCoNum+VSMALL));
    }
    else
    {
        CoNum = convectiveCoNum;
        maxCo = min(maxConvectiveCo,convectiveCoNum/(restrictInterfaceCoNum+VSMALL));
    }

    Info << "[DEBUG] CoNum: " << CoNum << endl;
    Info << "[DEBUG] maxCo: " << maxCo << endl;
