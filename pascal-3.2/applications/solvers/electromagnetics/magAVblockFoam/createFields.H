    // Magnetic properties
    Info << "Create magnetic properties" << nl << endl;
    IOdictionary magneticProperties
    (
        IOobject
        (
            "magneticProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar frequency (magneticProperties.lookup("frequency"));


    // Source current density (real part)
    Info << "Create field jsRe" << endl;
    volVectorField jsRe
    (
        IOobject
        (
            "jsRe",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Source current density (imaginary part)
    Info << "Create field jsIm" << nl << endl;
    volVectorField jsIm
    (
        IOobject
        (
            "jsIm",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    // Conductivity
    Info << "Create field sigma" << nl << endl;
    volScalarField sigma
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    // Magnetic vector potential (real part)
    Info << "Create field ARe" << nl << endl;
    volVectorField ARe
    (
        IOobject
        (
            "ARe",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Magnetic vector potential (imaginary part)
    Info << "Create field AIm" << nl << endl;
    volVectorField AIm
    (
        IOobject
        (
            "AIm",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    // Gauge potential (real part)
    Info << "Create field GRe" << endl;
    volScalarField GRe
    (
        IOobject
        (
            "GRe",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Gauge potential (imaginary part)
    Info << "Create field GIm" << nl << endl;
    volScalarField GIm
    (
        IOobject
        (
            "GIm",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Electric scalar potential (real part)
    Info << "Create field VRe" << endl;
    volScalarField VRe
    (
        IOobject
        (
            "VRe",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Electric scalar potential (imaginary part)
    Info << "Create field VIm " << nl << endl;
    volScalarField VIm
    (
        IOobject
        (
            "VIm",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    // Magnetic vector potential (complex)
    Info << "Create field AV" << nl << endl;
    volVector10Field AV
    (
        IOobject
        (
            "AV",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector10(word(),dimless,vector10::zero)
    );


    // Magnetic field density (real part)
    Info << "Create field BRe" << endl;
    volVectorField BRe
    (
        IOobject
        (
            "BRe",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::curl(ARe)
    );

    // Magnetic field density (imaginary part)
    Info << "Create field BIm" << nl << endl;
    volVectorField BIm
    (
        IOobject
        (
            "BIm",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::curl(AIm)
    );


    // Current density (real part)
    Info << "Create field jRe" << endl;
    volVectorField jRe
    (
        IOobject
        (
            "jRe",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector(word(),dimensionSet(0,-2,0,0,0,1,0),vector::zero)
    );

    // Current density (imaginary part)
    Info << "Create field jIm" << nl << endl;
    volVectorField jIm
    (
        IOobject
        (
            "jIm",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector(word(),dimensionSet(0,-2,0,0,0,1,0),vector::zero)
    );


    // Lorentz force
    Info << "Create field FL" << endl;
    volVectorField FL
    (
        IOobject
        (
            "FL",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensioned<vector>
        (
            "FL",
            dimMass/pow(dimLength,2)/pow(dimTime,2),
            pTraits<vector>::zero
        ),
        calculatedFvPatchField<vector>::typeName
    );
    Info << endl;


    // Lorentz force
    Info << "Create region field pB" << endl;
    volScalarField pB
    (
        IOobject
        (
            "pB",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensioned<scalar>
        (
            "pB",
            dimMass/dimLength/pow(dimTime,2),
            pTraits<scalar>::zero
        ),
        calculatedFvPatchField<scalar>::typeName
    );
    Info << endl;
