
    // Make fluxes relative
    fvc::makeRelative(phi,U);


    // Courant No
#   include "CourantNo.H"


    // Mesh Courant No
#   include "meshCourantNo.H"


    // Convective, combined Courant No
    scalar convectiveCoNum = max(CoNum, meshCoNum);


    // Interface CoNum
    scalar interfaceCoNum = interface.maxCourantNumber();

    Info << "Courant Number of tracked surface max: " << interfaceCoNum
        << endl;


    // Effective Courant No
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


    // Make fluxes absolute
    fvc::makeAbsolute(phi,U);


//    Info << "[DEBUG]: CoNum = " << CoNum << endl;
//    Info << "[DEBUG]: maxCo = " << maxCo << endl;
