
    // Region control
    Info << "Create region control" << endl;
    regionFvControl control_(mesh_);


    // Eddy current control
    Info << "Create eddy current control" << endl;
    eddyCurrentControl eddy(control_);

    Info << endl;
