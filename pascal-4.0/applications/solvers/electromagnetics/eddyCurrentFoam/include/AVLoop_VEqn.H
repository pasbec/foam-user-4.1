
    // TODO: Remove interfaceLabel? Use alias of calculatedFvPatchField or special patchField

    // Solve V
    if (control.meshIs3D())
    {
        SM_REGIONSCOPE(CONDUCTOR);

        volScalarField& sigma = globalStorage.sigma()[CONDUCTOR];

        label complexSign = 0;
        volScalarField* VPtr = NULL;
        volVectorField* VGradPtr = NULL;
        volVectorField* APtr = NULL;
        volVectorField* A0Ptr = NULL;
        surfaceScalarField* ddtAPhiPtr = NULL;

        if (control.switchV())
        {
            complexSign = 1;
            VPtr = &storage.VRe();
            VGradPtr = &globalStorage.VReGrad()[CONDUCTOR];
            APtr = &globalStorage.AIm()[CONDUCTOR];
            ddtAPhiPtr = &storage.ddtAImPhi();

            if (globalSettings.biotSavart)
            {
                A0Ptr = &globalStorage.A0Im()[CONDUCTOR];
            }
        }
        else
        {
            complexSign = -1;
            VPtr = &storage.VIm();
            VGradPtr = &globalStorage.VImGrad()[CONDUCTOR];
            APtr = &globalStorage.ARe()[CONDUCTOR];
            ddtAPhiPtr = &storage.ddtARePhi();

            if (globalSettings.biotSavart)
            {
                A0Ptr = &globalStorage.A0Re()[CONDUCTOR];
            }
        }

        volScalarField& V = *VPtr;
        volVectorField& VGrad = *VGradPtr;
        volVectorField& A = *APtr;
        surfaceScalarField& ddtAPhi = *ddtAPhiPtr;

// TODO: Use dedicated boundary condition?
        // Update gradient boundary condition for V
        {
            label patchI = control.interfaceLabel();

            if (isA<fixedGradientFvPatchScalarField>(V.boundaryField()[patchI]))
            {
                fixedGradientFvPatchScalarField& VP =
                    refCast<fixedGradientFvPatchScalarField>
                    (
                        V.boundaryField()[patchI]
                    );

                // Set correct gradient
                VP.gradient() = -ddtAPhi.boundaryField()[patchI]/VP.patch().magSf();
            }

            V.correctBoundaryConditions();
        }

        surfaceScalarField sigmaf("sigmaf", fvc::interpolate(sigma));

        // Read Laplacian scheme for discretising laplacian(sigmaf,V)
        word laplacianSchemeV
        (
            mesh.schemesDict().laplacianScheme
            (
                "laplacian(" + sigmaf.name() + "," + V.name() + ")"
            )
        );

        // Prepare V-equation
        fvScalarMatrix VEqn(fvm::laplacian(sigmaf, V));

        if (laplacianSchemeV != "jumpGauss")
        {
            // NOTE: We do NOT use sigmaf here!
            volVectorField sigmaGrad = fvc::grad(sigma);

            // Explicit source term for V due to variable sigma
            // in conducting region with implied condition div(A)=0
            VEqn -= complexSign * omega0 * (sigmaGrad & A);
//          VEqn -= complexSign * omega0 * sigma * fvc::div(A);

            if (globalSettings.biotSavart)
            {
                volVectorField& A0 = *A0Ptr;

                VEqn -= complexSign * omega0 * (sigmaGrad & A0);
            }
        }

        // Read solver for V
        word solverV
        (
            mesh.solutionDict().solver(V.name()).lookup("solver")
        );

        // Solve V-equation
        if(solverV == "deflatedPCG")
        {
            // Currently a deflated solver is necessary.
            // Therefor the matrix needs to be kept singluar!
            VEqn.solve(control.subDict(Region::CONDUCTOR, V.name()));
        }
        else
        {
            FatalErrorIn
            (
                Control::typeName + " | VEqn.H :\n"
            )   << "Only the 'deflatedPCG' solver is currently supported for"
                << " the solution of VRe/VIm."
                << abort(FatalError);
        }

        // Correct boundary conditions
        V.correctBoundaryConditions();

// TODO: Use relaxation?
//         VReGrad.storePrevIter();
//         VImGrad.storePrevIter();

        if (laplacianSchemeV != "jumpGauss")
        {
            VGrad = fvc::reconstruct(VEqn.flux()/sigmaf);
        }
        else
        {
            VGrad = fvc::grad(V);
        }

// TODO: Use relaxation?
//         VReGrad.relax(0.5);
//         VImGrad.relax(0.5);
    }
