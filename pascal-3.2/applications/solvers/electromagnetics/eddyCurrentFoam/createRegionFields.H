
# define setRegionReference(regionI, vf)                                      \
                                                                              \
    typeof(vf##_[regionI])& vf = vf##_[regionI]; if (&vf != &vf) {}

# define setRegion(regionI)                                                   \
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
