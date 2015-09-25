
    // Assemble matrix for A
#   include "AEqn.H"

    // Assemble matrix for V
#   include "VEqn.H"

    // Solve AV-system
    if (AVmethod == "segregated")
    {
        int nSegregatedCorr = AVmethodDict.lookupOrDefault<int>("corr",10);

        for
        (
            int segregatedCorr=0;
            segregatedCorr<nSegregatedCorr;
            segregatedCorr++
        )
        {
#           include "solveAVsegregated.H"
        }
    }
    else if (AVmethod == "coupled")
    {
        int nCoupledCorr = AVmethodDict.lookupOrDefault<int>("corr",4);

        for
        (
            int coupledCorr=0;
            coupledCorr<nCoupledCorr;
            coupledCorr++
        )
        {
#           include "solveAVcoupled.H"
        }
    }
    else
    {
        FatalErrorIn("AVEqn.H")
        << "Method for solving AV-system '" << AVmethod
        << "' is not supported"
        << abort(FatalError);
    }

    // Update corrected gradient of V
    {
        // Interpolate V to face centers
        surfaceScalarField sVRe
        (
            "sVRe",
            fvc::interpolate(VRe)
        );
        surfaceScalarField sVIm
        (
            "sVIm",
            fvc::interpolate(VIm)
        );

        // Enforce zero normal gradient for
        // V (d(V)/dn=0)
        {
            // Mesh cell-face addressing
            const unallocLabelList& own = mesh.owner();
            const unallocLabelList& nei = mesh.neighbour();

            // Replace face values with cell values
            // from conducting side
            forAll (cnPNfaces, facei)
            {
                label faceI = cnPNfaces[facei];
                sVRe[faceI] = VRe[own[faceI]];
                sVIm[faceI] = VIm[own[faceI]];
            }
            forAll (cnNPfaces, facei)
            {
                label faceI = cnNPfaces[facei];
                sVRe[faceI] = VRe[nei[faceI]];
                sVIm[faceI] = VIm[nei[faceI]];
            }
        }

        // Update corrected gradient of V
        VReGrad = fvc::grad(sVRe);
        VImGrad = fvc::grad(sVIm);

        // Remove all gradient values from
        // non-conducting region
        forAll(nCells, celli)
        {
            label cellI = nCells[celli];
            VReGrad[cellI] = vector::zero;
            VImGrad[cellI] = vector::zero;
        }
    }



    // Calculate magnetic field density
    BRe == fvc::curl(ARe);
    BIm == fvc::curl(AIm);

    // Calculate eddy current density
    jRe ==   alpha * AIm - sigma * VReGrad;
    jIm == - alpha * ARe - sigma * VImGrad;



    // Gauge fixing: Solve for gauge fixing potential G
    {
        // Assemble G-equations
        fvScalarMatrix GReEqn
        (
            fvm::laplacian(GRe)
          + fvc::div(ARe)
        );
        fvScalarMatrix GImEqn
        (
            fvm::laplacian(GIm)
          + fvc::div(AIm)
        );

        // Reference settings
        label refCell = mesh.findCell(point(0.00, 0.00, 0.01));
        scalar refValue = 0.0;

        GReEqn.setReference(refCell, refValue);
        GImEqn.setReference(refCell, refValue);

        // Solve G
        GReEqn.solve();
        GImEqn.solve();

        GRe.correctBoundaryConditions();
        GIm.correctBoundaryConditions();

        // TODO: Remove after debug
        volVectorField AReCoulomb("AReCoulomb", ARe + fvc::grad(GRe));
        volVectorField AImCoulomb("AImCoulomb", AIm + fvc::grad(GIm));
        AReCoulomb.write();
        AImCoulomb.write();
    }

