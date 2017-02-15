
    // Solve A0 in conducting region
    if (globalSettings.biotSavart)
    {
        SM_REGIONSCOPE(CONDUCTOR);

        volScalarField& nur = globalStorage.nur()[CONDUCTOR];

        volVectorField& A0Re = globalStorage.A0Re()[CONDUCTOR];
        volVectorField& A0Im = globalStorage.A0Im()[CONDUCTOR];

        // Interpolate relative magnetic reluctivity
        surfaceScalarField nurf("nurf", fvc::interpolate(nur));

        // Assemble base-matrices for A0
        fvVectorMatrix A0ReEqn
        (
          - fvm::laplacian(nurf,A0Re)
        );

        fvVectorMatrix A0ImEqn
        (
          - fvm::laplacian(nurf,A0Im)
        );

        solve(A0ReEqn);
        solve(A0ImEqn);
    }
