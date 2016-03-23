
# define setRegionScope(regionI)                                              \
                                                                              \
    setConstRegionReference(regionI, fvMesh, mesh);                           \
    setRegionReference(regionI, volScalarField, sigma);                       \
    setRegionReference(regionI, volVectorField, ARe);                         \
    setRegionReference(regionI, volVectorField, AIm);                         \
    setRegionReference(regionI, volVectorField, VReGrad);                     \
    setRegionReference(regionI, volVectorField, VImGrad);                     \
    setRegionReference(regionI, volVectorField, FL);                          \
    setRegionReference(regionI, volScalarField, pB);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Control class
    Info << "Create control" << nl << endl;
    eddyCurrentControl control(args, runTime, mesh_);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Conductivity
    infoCreateRegionField(sigma);
    regionVolScalarField sigma_
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    );
    Info << nl;


    // Magnetic vector potential (real part)
    infoCreateRegionField(ARe);
    regionVolVectorField ARe_
    (
        IOobject
        (
            "ARe",
            runTime.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    );

    // Magnetic vector potential (imaginary part)
    infoCreateRegionField(AIm);
    regionVolVectorField AIm_
    (
        IOobject
        (
            "AIm",
            runTime.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    );
    Info << nl;


    // Gradient of electric scalar potential (real part)
    infoCreateRegionField(VReGrad);
    regionVolVectorField VReGrad_
    (
        IOobject
        (
            "VReGrad",
            runTime.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensioned<vector>
        (
            word(),
            dimensionSet(1,1,-3,0,0,-1,0),
            pTraits<vector>::zero
        ),
        calculatedFvPatchField<vector>::typeName
    );

    // Gradient of electric scalar potential (imaginary part)
    infoCreateRegionField(VImGrad);
    regionVolVectorField VImGrad_
    (
        IOobject
        (
            "VImGrad",
            runTime.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensioned<vector>
        (
            word(),
            dimensionSet(1,1,-3,0,0,-1,0),
            pTraits<vector>::zero
        ),
        calculatedFvPatchField<vector>::typeName
    );
    Info << nl;


    // Lorentz force
    infoCreateRegionField(FL);
    regionVolVectorField FL_
    (
        IOobject
        (
            "FL",
            runTime.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensioned<vector>
        (
            word(),
            dimMass/pow(dimLength,2)/pow(dimTime,2),
            pTraits<vector>::zero
        ),
        calculatedFvPatchField<vector>::typeName
    );
    Info << endl;


    // Lorentz force
    infoCreateRegionField(pB);
    regionVolScalarField pB_
    (
        IOobject
        (
            "pB",
            runTime.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensioned<scalar>
        (
            word(),
            dimMass/dimLength/pow(dimTime,2),
            pTraits<scalar>::zero
        ),
        calculatedFvPatchField<scalar>::typeName
    );
    Info << endl;
