
    // Electric scalar potential (real part)
    infoCreateField(eddy.conductor(), VRe);
    volScalarField VRe
    (
        IOobject
        (
            "VRe",
            runTime.timeName(),
            mesh_[eddy.conductor()],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_[eddy.conductor()]
    );
    uniformDimensionedScalarField VReAverage
    (
        IOobject
        (
            "VReAverage",
            runTime.timeName(),
            mesh_[eddy.conductor()],
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
    infoCreateField(eddy.conductor(), VIm);
    volScalarField VIm
    (
        IOobject
        (
            "VIm",
            runTime.timeName(),
            mesh_[eddy.conductor()],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_[eddy.conductor()]
    );
    uniformDimensionedScalarField VImAverage
    (
        IOobject
        (
            "VImAverage",
            runTime.timeName(),
            mesh_[eddy.conductor()],
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        dimensionedScalar(
            word(),
            VIm.dimensions(),
            gAverage(VIm.internalField())
        )
    );
    Info << nl;
