
    // Electric scalar potential (real part)
    Info << "Create region field VRe" << endl;
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

    // Electric scalar potential (imaginary part)
    Info << "Create region field VIm" << endl;
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


    // Gradient of electric scalar potential (real part)
    Info << "Create region field VReGrad" << endl;
    regionVolVectorField VReGrad_
    (
        IOobject
        (
            "VReGrad",
            runTime.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(VRe_)
    );

    // Gradient of electric scalar potential (imaginary part)
    Info << "Create region field VImGrad" << endl;
    regionVolVectorField VImGrad_
    (
        IOobject
        (
            "VImGrad",
            runTime.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(VIm_)
    );
    Info << endl;
