
# define setRegionReference(region, refName)                        \
                                                                    \
    typeof(refName##_[region])& refName =                           \
        refName##_[region];                                         \
                                                                    \
    if (&refName != &refName) {} // Hack to trick -Wunused-variable

# define setRegionReferences(region)                                \
                                                                    \
    setRegionReference(region, mesh);                               \
    setRegionReference(region, sigma);                              \
    setRegionReference(region, ARe);                                \
    setRegionReference(region, AIm);                                \
    setRegionReference(region, VReGrad);                            \
    setRegionReference(region, VImGrad);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Conductivity
    Info << "Create region field sigma" << nl << endl;
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


    // Magnetic vector potential (real part)
    Info << "Create region field ARe" << endl;
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
    Info << "Create region field AIm" << nl << endl;
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
    Info << "Create region field VImGrad" << nl << endl;
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
        mesh_,
        dimensioned<vector>
        (
            word(),
            dimensionSet(1,1,-3,0,0,-1,0),
            pTraits<vector>::zero
        ),
        calculatedFvPatchField<vector>::typeName
    );
