

    // Source current density (real part)
    infoCreateField(eddy.base(), jsRe);
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
    infoCreateField(eddy.base(), jsIm);
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
    Info << nl;

    // Magnetic field density (real part)
    infoCreateField(eddy.base(), BRe);
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
    infoCreateField(eddy.base(), BIm);
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
    Info << nl;


    // Current density (real part)
    infoCreateField(eddy.base(), jRe);
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
    infoCreateField(eddy.base(), jIm);
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
    Info << nl;


    // Lorentz force
    infoCreateField(eddy.base(), FL);
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
    Info << nl;


    // Lorentz force
    infoCreateField(eddy.base(), pB);
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
    Info << nl;
