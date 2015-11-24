

    // Source current density (real part)
    Info << "Create field jsRe" << endl;
    volVectorField jsRe
    (
        IOobject
        (
            "jsRe",
            runTime.timeName(),
            mesh_[eddy.base()],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_[eddy.base()]
    );

    // Source current density (imaginary part)
    Info << "Create field jsIm" << nl << endl;
    volVectorField jsIm
    (
        IOobject
        (
            "jsIm",
            runTime.timeName(),
            mesh_[eddy.base()],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_[eddy.base()]
    );

    // Magnetic field density (real part)
    Info << "Create field BRe" << endl;
    volVectorField BRe
    (
        IOobject
        (
            "BRe",
            runTime.timeName(),
            mesh_[eddy.base()],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_[eddy.base()]
    );

    // Magnetic field density (imaginary part)
    Info << "Create field BIm" << nl << endl;
    volVectorField BIm
    (
        IOobject
        (
            "BIm",
            runTime.timeName(),
            mesh_[eddy.base()],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_[eddy.base()]
    );


    // Current density (real part)
    Info << "Create field jRe" << endl;
    volVectorField jRe
    (
        IOobject
        (
            "jRe",
            runTime.timeName(),
            mesh_[eddy.base()],
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_[eddy.base()],
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
            mesh_[eddy.base()],
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_[eddy.base()],
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
            mesh_[eddy.base()],
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_[eddy.base()],
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
    Info << "Create field pB" << endl;
    volScalarField pB
    (
        IOobject
        (
            "pB",
            runTime.timeName(),
            mesh_[eddy.base()],
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_[eddy.base()],
        dimensioned<scalar>
        (
            "pB",
            dimMass/dimLength/pow(dimTime,2),
            pTraits<scalar>::zero
        ),
        calculatedFvPatchField<scalar>::typeName
    );
    Info << endl;
