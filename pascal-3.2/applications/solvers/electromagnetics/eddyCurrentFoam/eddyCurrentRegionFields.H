
# define setRegionScope(regionI)                                              \
                                                                              \
    setRegionReference(regionI, fvMesh, mesh);                                \
    setRegionReference(regionI, volScalarField, sigma);                       \
    setRegionReference(regionI, volVectorField, ARe);                         \
    setRegionReference(regionI, volVectorField, AIm);                         \
    setRegionReference(regionI, volVectorField, VReGrad);                     \
    setRegionReference(regionI, volVectorField, VImGrad);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Control class
    Info << "Create control" << nl << endl;
    eddyCurrentControl control(mesh_);

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
