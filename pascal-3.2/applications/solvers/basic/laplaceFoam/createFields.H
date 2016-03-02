    Info<< "Reading field psi\n" << endl;

    volScalarField psi
    (
        IOobject
        (
            "psi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Creating field psiGrad\n" << endl;

    volVectorField psiGrad
    (
        IOobject
        (
            "psiGrad",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(psi)
    );


    Info<< "Creating field psiSnGradGrad\n" << endl;

    volVectorField psiSnGradGrad
    (
        IOobject
        (
            "psiSnGradGrad",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::reconstruct(fvc::snGrad(psi) * mesh.magSf())
    );


    Info<< "Creating field psiPhi\n" << endl;

    surfaceScalarField psiPhi
    (
        IOobject
        (
            "psiPhi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            word(),
            psi.dimensions()*dimLength,
            0
        ),
        calculatedFvPatchVectorField::typeName
    );


    Info<< "Creating field psiFlux\n" << endl;

    volVectorField psiFlux
    (
        IOobject
        (
            "psiFlux",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            word(),
            psi.dimensions()/dimLength,
            vector::zero
        ),
        calculatedFvPatchVectorField::typeName
    );


    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    Info<< "Reading diffusivity psiGamma\n" << endl;

    dimensionedScalar psiGamma
    (
        transportProperties.lookup("psiGamma")
    );

    label psiRefCell = 0;
    scalar psiRefValue = 0.0;
    setRefCell(psi, mesh.solutionDict().subDict("SIMPLE"), psiRefCell, psiRefValue);
