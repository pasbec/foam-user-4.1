

    // Source current density (real part)
    Info << "Create field jsRe" << endl;
    volVectorField jsRe
    (
        IOobject
        (
            "jsRe",
            runTime.timeName(),
            mesh_[eddy.baseRegion()],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_[eddy.baseRegion()]
    );

    // Source current density (imaginary part)
    Info << "Create field jsIm" << nl << endl;
    volVectorField jsIm
    (
        IOobject
        (
            "jsIm",
            runTime.timeName(),
            mesh_[eddy.baseRegion()],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_[eddy.baseRegion()]
    );


//     // Electric scalar potential (real part)
//     Info << "Create field VRe" << endl;
//     volScalarField VRe
//     (
//         IOobject
//         (
//             "VRe",
//             runTime.timeName(),
//             mesh_[eddy.baseRegion()],
//             IOobject::MUST_READ,
//             IOobject::AUTO_WRITE
//         ),
//         mesh_[eddy.baseRegion()]
//     );
//     uniformDimensionedScalarField VReAverage
//     (
//         IOobject
//         (
//             "VReAverage",
//             runTime.timeName(),
//             mesh_[eddy.baseRegion()],
//             IOobject::READ_IF_PRESENT,
//             IOobject::AUTO_WRITE
//         ),
//         dimensionedScalar(
//             word(),
//             VRe.dimensions(),
//             gAverage(VRe.internalField())
//         )
//     );

//     // Electric scalar potential (imaginary part)
//     Info << "Create field VIm " << nl << endl;
//     volScalarField VIm
//     (
//         IOobject
//         (
//             "VIm",
//             runTime.timeName(),
//             mesh_[eddy.baseRegion()],
//             IOobject::MUST_READ,
//             IOobject::AUTO_WRITE
//         ),
//         mesh_[eddy.baseRegion()]
//     );
//     uniformDimensionedScalarField VImAverage
//     (
//         IOobject
//         (
//             "VImAverage",
//             runTime.timeName(),
//             mesh_[eddy.baseRegion()],
//             IOobject::READ_IF_PRESENT,
//             IOobject::AUTO_WRITE
//         ),
//         dimensionedScalar(
//             word(),
//             VIm.dimensions(),
//             gAverage(VIm.internalField())
//         )
//     );


    // Magnetic field density (real part)
    Info << "Create field BRe" << endl;
    volVectorField BRe
    (
        IOobject
        (
            "BRe",
            runTime.timeName(),
            mesh_[eddy.baseRegion()],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_[eddy.baseRegion()]
    );

    // Magnetic field density (imaginary part)
    Info << "Create field BIm" << nl << endl;
    volVectorField BIm
    (
        IOobject
        (
            "BIm",
            runTime.timeName(),
            mesh_[eddy.baseRegion()],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_[eddy.baseRegion()]
    );


    // Current density (real part)
    Info << "Create field jRe" << endl;
    volVectorField jRe
    (
        IOobject
        (
            "jRe",
            runTime.timeName(),
            mesh_[eddy.baseRegion()],
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_[eddy.baseRegion()],
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
            mesh_[eddy.baseRegion()],
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_[eddy.baseRegion()],
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
            mesh_[eddy.baseRegion()],
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_[eddy.baseRegion()],
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
            mesh_[eddy.baseRegion()],
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_[eddy.baseRegion()],
        dimensioned<scalar>
        (
            "pB",
            dimMass/dimLength/pow(dimTime,2),
            pTraits<scalar>::zero
        ),
        calculatedFvPatchField<scalar>::typeName
    );
    Info << endl;
