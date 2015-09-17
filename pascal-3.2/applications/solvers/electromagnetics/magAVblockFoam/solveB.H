
    // Calculate msigma
    volScalarField msigma
    (
        "msigma",
        sigma
    );

    // Manipulate msigma:
    //
    //  We will replace "zero"-values with very small
    //  values at all faces in order to avoid a singular
    //  solution
    {
        // Parameters
        scalar msigmaMin = 1.0;

        // Manipulate internal faces
        forAll(msigma.internalField(), cellI)
        {
            if
            (msigma.internalField()[cellI] < msigmaMin)
            {
                msigma.internalField()[cellI] = msigmaMin;
            }
        }

        // Manipulate boundary faces
        msigma.boundaryField() == msigmaMin;
    }

    // Interpolate sigma to face centers
    surfaceScalarField msigmaf
    (
        "msigmaf",
        fvc::interpolate(msigma, "interpolate(sigma)")
    );

    // Create surface field for magnetic diffusivity
    surfaceScalarField mgammaf
    (
        "gammaf",
        rMu0/msigmaf
    );

    {
        // Real part
        fvVectorMatrix BReEqn
        (
            mu0 * fvm::laplacian(mgammaf, BRe, "laplacian(BRe)")
          + mu0 * fvc::laplacian(mgammaf, fvc::grad(GRe), "laplacian(BRe)")
//           + mu0 * omega * BIm
          + mu0 * omega * fvc::grad(GIm)
          + mu0 * fvc::curl(jsRe/msigma)
        );

        // Imaginary part
        fvVectorMatrix BImEqn
        (
            mu0 * fvm::laplacian(mgammaf, BIm, "laplacian(BIm)")
          + mu0 * fvc::laplacian(mgammaf, fvc::grad(GIm), "laplacian(BIm)")
//           - mu0 * omega * BRe
          - mu0 * omega * fvc::grad(GRe)
          + mu0 * fvc::curl(jsIm/msigma)
        );

        // Prepare block system
        fvBlockMatrix<vector6> BEqn(B);

        // Insert equations for B
        BEqn.insertEquation(0, BReEqn);
        BEqn.insertEquation(3, BImEqn);

        // Coupling coeffs
        scalarField alpha(mesh.nCells(), mu0.value() * omega.value());

        // Add coupling terms for BRe-BIm coupling
        BEqn.insertEquationCoupling(0, 3,  alpha);
        BEqn.insertEquationCoupling(3, 0, -alpha);
        BEqn.insertEquationCoupling(1, 4,  alpha);
        BEqn.insertEquationCoupling(4, 1, -alpha);
        BEqn.insertEquationCoupling(2, 5,  alpha);
        BEqn.insertEquationCoupling(5, 2, -alpha);

        // Block coupled solver call to solve B
        BEqn.solve();

        // Retrieve solution
        BEqn.retrieveSolution(0, BRe.internalField());
        BEqn.retrieveSolution(3, BIm.internalField());

        BRe.correctBoundaryConditions();
        BIm.correctBoundaryConditions();

        fvScalarMatrix GReEqn
        (
            fvm::laplacian(GRe)
          - fvc::div(BRe)
        );

        fvScalarMatrix GImEqn
        (
            fvm::laplacian(GIm)
          - fvc::div(BIm)
        );

        GReEqn.solve();
        GImEqn.solve();

        GRe.correctBoundaryConditions();
        GIm.correctBoundaryConditions();

//         BRe += fvc::grad(GRe);
//         BIm += fvc::grad(GIm);
// //
//         BRe.correctBoundaryConditions();
//         BIm.correctBoundaryConditions();
    }

//     jRe == rMu0 * fvc::curl(BRe) - jsRe;
//     jIm == rMu0 * fvc::curl(BIm) - jsIm;

    jRe == rMu0 * fvc::curl(BRe);
    jIm == rMu0 * fvc::curl(BIm);

    volScalarField BReDiv("BReDiv", fvc::div(BRe)); BReDiv.write();
    volScalarField BImDiv("BImDiv", fvc::div(BIm)); BImDiv.write();
