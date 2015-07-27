    // Magnetic properties
    Info << "Create magnetic properties for"
        << " region " << mesh_[defaultRegionID].name()
        << nl << endl;
    IOdictionary magneticProperties
    (
        IOobject
        (
            "magneticProperties",
            runTime.constant(),
            mesh_[defaultRegionID],
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar frequency (magneticProperties.lookup("frequency"));


    // Source current density (real part)
    Info << "Create field jsRe for"
        << " region " << mesh_[defaultRegionID].name()
        << endl;
    volVectorField jsRe
    (
        IOobject
        (
            "jsRe",
            runTime.timeName(),
            mesh_[defaultRegionID],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_[defaultRegionID]
    );

    // Source current density (imaginary part)
    Info << "Create field jsIm for"
        << " region " << mesh_[defaultRegionID].name()
        << nl << endl;
    volVectorField jsIm
    (
        IOobject
        (
            "jsIm",
            runTime.timeName(),
            mesh_[defaultRegionID],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_[defaultRegionID]
    );


    // Conductivity
    Info << "Create field sigma for"
        << " region " << mesh_[defaultRegionID].name()
        << nl << endl;
    volScalarField sigma
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            mesh_[defaultRegionID],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_[defaultRegionID]
    );


    // Megnetic vector potential (real part)
    Info << "Create field ARe for"
        << " region " << mesh_[defaultRegionID].name()
        << endl;
    volVectorField ARe
    (
        IOobject
        (
            "ARe",
            runTime.timeName(),
            mesh_[defaultRegionID],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_[defaultRegionID]
    );

    // Megnetic vector potential (imaginary part)
    Info << "Create field AIm for"
        << " region " << mesh_[defaultRegionID].name()
        << endl;
    volVectorField AIm
    (
        IOobject
        (
            "AIm",
            runTime.timeName(),
            mesh_[defaultRegionID],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_[defaultRegionID]
    );

    // Megnetic vector potential (complex)
    Info << "Create field A for"
        << " region " << mesh_[defaultRegionID].name()
        << nl << endl;
    volVector6Field A
    (
        IOobject
        (
            "A",
            runTime.timeName(),
            mesh_[defaultRegionID],
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_[defaultRegionID],
        dimensionedVector6(word(),dimless,vector6::zero)
    );


    // Magnetic field density (real part)
    Info << "Create field BRe for"
        << " region " << mesh_[defaultRegionID].name()
        << endl;
    volVectorField BRe
    (
        IOobject
        (
            "BRe",
            runTime.timeName(),
            mesh_[defaultRegionID],
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::curl(ARe)
    );

    // Magnetic field density (imaginary part)
    Info << "Create field BIm for"
        << " region " << mesh_[defaultRegionID].name()
        << nl << endl;
    volVectorField BIm
    (
        IOobject
        (
            "BIm",
            runTime.timeName(),
            mesh_[defaultRegionID],
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::curl(AIm)
    );


    // Current density (real part)
    Info << "Create field jRe for"
        << " region " << mesh_[defaultRegionID].name()
        << endl;
    volVectorField jRe
    (
        IOobject
        (
            "jRe",
            runTime.timeName(),
            mesh_[defaultRegionID],
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_[defaultRegionID],
        dimensionedVector(word(),dimensionSet(0,-2,0,0,0,1,0),vector::zero)
    );

    // Current density (imaginary part)
    Info << "Create field jIm for"
        << " region " << mesh_[defaultRegionID].name()
        << nl << endl;
    volVectorField jIm
    (
        IOobject
        (
            "jIm",
            runTime.timeName(),
            mesh_[defaultRegionID],
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_[defaultRegionID],
        dimensionedVector(word(),dimensionSet(0,-2,0,0,0,1,0),vector::zero)
    );
