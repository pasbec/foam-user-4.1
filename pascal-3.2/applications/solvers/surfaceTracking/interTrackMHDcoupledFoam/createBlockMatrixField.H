
    // FIXME [High]: Why can't we use symmetryPlane as patch type???
    Info << "Creating field Up" << nl << endl;
    volVector4Field Up
    (
        IOobject
        (
        "Up",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector4("zero", dimless, vector4::zero)
    );