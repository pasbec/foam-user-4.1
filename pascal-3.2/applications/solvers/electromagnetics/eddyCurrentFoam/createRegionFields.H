
# define infoCreateRegionField(f)                                             \
                                                                              \
    Info << "Create region field " << #f << endl;

# define infoCreateField(regionI, f)                                          \
                                                                              \
    Info << "Create field " << #f                                             \
        << " (" << mesh_.regionName(regionI) << ")" << endl;

# define setRegionReference(regionI, f)                                       \
                                                                              \
    typeof(f##_[regionI])& f = f##_[regionI]; if (&f != &f) {}

# define setRegionScope(regionI)                                              \
                                                                              \
    setRegionReference(regionI, mesh);                                        \
    setRegionReference(regionI, sigma);                                       \
    setRegionReference(regionI, ARe);                                         \
    setRegionReference(regionI, AIm);                                         \
    setRegionReference(regionI, VReGrad);                                     \
    setRegionReference(regionI, VImGrad);

# define updateSigmaA(regionI)                                                \
                                                                              \
    sigma_.mapExtrapolate(regionI);                                           \
    ARe_.mapInterpolate(regionI);                                             \
    AIm_.mapInterpolate(regionI);

# define updateRmapGradient(regionI,vsf)                                      \
                                                                              \
    vsf##Grad_[regionI] = fvc::grad(vsf);                                     \
    vsf##Grad_.rmap(regionI);

# define updateGradientV(regionI)                                             \
                                                                              \
    updateRmapGradient(regionI,VRe);                                          \
    updateRmapGradient(regionI,VIm);


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
    infoCreateRegionField(VImGrad);
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
    Info << nl;
