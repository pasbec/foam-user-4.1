
    // Electric scalar potential (real part)
    infoCreateField(control.conductor(), VRe);
    volScalarField VRe
    (
        IOobject
        (
            "VRe",
            runTime.timeName(),
            mesh_[control.conductor()],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_[control.conductor()]
    );
    uniformDimensionedScalarField VReAverage
    (
        IOobject
        (
            "VReAverage",
            runTime.timeName(),
            mesh_[control.conductor()],
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
    infoCreateField(control.conductor(), VIm);
    volScalarField VIm
    (
        IOobject
        (
            "VIm",
            runTime.timeName(),
            mesh_[control.conductor()],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_[control.conductor()]
    );
    uniformDimensionedScalarField VImAverage
    (
        IOobject
        (
            "VImAverage",
            runTime.timeName(),
            mesh_[control.conductor()],
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
