
    // Electric scalar potential (real part)
    Info << "Create region field VRe"
        << endl;
    regionVolScalarField VRe_
    (
        IOobject
        (
            "VRe",
            runTime.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    );
    Info << endl;

    // Electric scalar potential (real part)
    Info << "Create region field VIm"
        << endl;
    regionVolScalarField VIm_
    (
        IOobject
        (
            "VIm",
            runTime.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    );
    Info << endl;
