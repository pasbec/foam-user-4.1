
    // Electric scalar potential (real part)
    Info << "Create field VRe" << endl;
    volScalarField VRe
    (
        IOobject
        (
            "VRe",
            runTime.timeName(),
            mesh_[eddy.conductorRegion()],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_[eddy.conductorRegion()]
    );
    uniformDimensionedScalarField VReAverage
    (
        IOobject
        (
            "VReAverage",
            runTime.timeName(),
            mesh_[eddy.conductorRegion()],
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        dimensionedScalar(
            word(),
            VRe.dimensions(),
            gAverage(VRe.internalField())
        )
    );

    // Electric scalar potential (imaginary part)
    Info << "Create field VIm " << nl << endl;
    volScalarField VIm
    (
        IOobject
        (
            "VIm",
            runTime.timeName(),
            mesh_[eddy.conductorRegion()],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_[eddy.conductorRegion()]
    );
    uniformDimensionedScalarField VImAverage
    (
        IOobject
        (
            "VImAverage",
            runTime.timeName(),
            mesh_[eddy.conductorRegion()],
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        dimensionedScalar(
            word(),
            VIm.dimensions(),
            gAverage(VIm.internalField())
        )
    );
