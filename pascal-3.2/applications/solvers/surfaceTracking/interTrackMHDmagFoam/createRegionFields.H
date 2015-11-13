
    // Lorentz force
    Info << "Create region field FL" << endl;
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
    Info << "Create region field pB" << endl;
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


//     // New points for mesh movement
//     Info << "Create region field FL"
//         << endl;
//     regionPointVectorField newPoints_
//     (
//         IOobject
//         (
//             "newPoints",
//             runTime.timeName(),
//             mesh_,
//             IOobject::READ_IF_PRESENT,
//             IOobject::AUTO_WRITE
//         ),
//         pointMesh_
//     );
//     Info << endl;
