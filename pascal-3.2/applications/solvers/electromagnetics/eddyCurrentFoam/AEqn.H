
    // Assemble base-matrices for A
    fvVectorMatrix AReEqn
    (
        rMu0 * fvm::laplacian(ARe)
//    + alpha * AIm
//    - sigma * VReGrad
      + jsRe
    );
    fvVectorMatrix AImEqn
    (
        rMu0 * fvm::laplacian(AIm)
//    - alpha * ARe
//    - sigma * VImGrad
      + jsIm
    );



    // Create list with boundary types
    wordList ARePatches(mesh.boundary().size());
    wordList AImPatches(mesh.boundary().size());

    forAll(ARe.boundaryField(), patchI)
    {
        ARePatches[patchI] =
            ARe.boundaryField()[patchI].type();
    }
    forAll(AIm.boundaryField(), patchI)
    {
        AImPatches[patchI] =
            AIm.boundaryField()[patchI].type();
    }

    // TODO: Check, that ARe/AIm use correct boundary
    //       conditions! Currently, only tangentialMagnetic
    //       is supported!

    // Create component fields for A
    volScalarField AxRe
    (
        IOobject
        (
            "AxRe",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        ),
        ARe.component(0) + dimensionedScalar(word(),ARe.dimensions(),VSMALL),
        ARePatches
    );
    volScalarField AyRe
    (
        IOobject
        (
            "AyRe",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        ),
        ARe.component(1) + dimensionedScalar(word(),ARe.dimensions(),VSMALL),
        ARePatches
    );
    volScalarField AzRe
    (
        IOobject
        (
            "AzRe",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        ),
        ARe.component(2) + dimensionedScalar(word(),ARe.dimensions(),VSMALL),
        ARePatches
    );
    volScalarField AxIm
    (
        IOobject
        (
            "AxIm",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        ),
        AIm.component(0) + dimensionedScalar(word(),AIm.dimensions(),VSMALL),
        AImPatches
    );
    volScalarField AyIm
    (
        IOobject
        (
            "AyIm",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        ),
        AIm.component(1) + dimensionedScalar(word(),AIm.dimensions(),VSMALL),
        AImPatches
    );
    volScalarField AzIm
    (
        IOobject
        (
            "AzIm",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        ),
        AIm.component(2) + dimensionedScalar(word(),AIm.dimensions(),VSMALL),
        AImPatches
    );



    // Assemble A-equations
    fvScalarMatrix AxReEqn
    (
        rMu0 * fvm::laplacian(AxRe)
//    + alpha * AxIm
//    - sigma * VReGrad.component(0)
      + jsRe.component(0)
    );
    fvScalarMatrix AyReEqn
    (
        rMu0 * fvm::laplacian(AyRe)
//    + alpha * AyIm
//    - sigma * VReGrad.component(1)
      + jsRe.component(1)
    );
    fvScalarMatrix AzReEqn
    (
        rMu0 * fvm::laplacian(AzRe)
//    + alpha * AzIm
//    - sigma * VReGrad.component(2)
      + jsRe.component(2)
    );
    fvScalarMatrix AxImEqn
    (
        rMu0 * fvm::laplacian(AxIm)
//    - alpha * AxRe
//    - sigma * VImGrad.component(0)
      + jsIm.component(0)
    );
    fvScalarMatrix AyImEqn
    (
        rMu0 * fvm::laplacian(AyIm)
//    - alpha * AyRe
//    - sigma * VImGrad.component(1)
      + jsIm.component(1)
    );
    fvScalarMatrix AzImEqn
    (
        rMu0 * fvm::laplacian(AzIm)
//    - alpha * AzRe
//    - sigma * VImGrad.component(2)
      + jsIm.component(2)
    );
