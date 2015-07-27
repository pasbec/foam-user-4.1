
    rAU = 1.0/UEqn.A();

    volVectorField*     HbyAPtr    = NULL;
    surfaceScalarField* rAUfPtr    = NULL;
    surfaceScalarField* phiFLPtr   = NULL;
    surfaceScalarField* phipBPtr   = NULL;
    surfaceScalarField* phiHbyAPtr = NULL;
    surfaceScalarField* phiRecPtr  = NULL;

    if (lorentzForceNonAlberto)
    {
        U = rAU * UEqn.H();

        phi = (fvc::interpolate(U) & mesh.Sf());
    }
    else
    {
        // Create HbyA from U to get BC from U
        HbyAPtr = new volVectorField
        (
            "HbyA",
            U
        );

        // Re-calculate inner field of HbyA
        *HbyAPtr = rAU * UEqn.H();

        rAUfPtr = new surfaceScalarField
        (
            "rAUf",
            fvc::interpolate(rAU)
        );

        phiFLPtr = new surfaceScalarField
        (
            "phiFL",
            *rAUfPtr * ( fvc::interpolate(rho*FL) & mesh.Sf() ) * lorentzForceVolumeFactor
        );

        phipBPtr = new surfaceScalarField
        (
            "phipB",
            *rAUfPtr * ( fvc::snGrad(pB) * mesh.magSf() ) * lorentzForceRotationalFactor * lorentzForceVolumeFactor
        );

        phiHbyAPtr = new surfaceScalarField
        (
            "phiHbyA",
            ( fvc::interpolate(*HbyAPtr) & mesh.Sf() )
          + *phiFLPtr
          + *phipBPtr
        );

        phiRecPtr = new surfaceScalarField
        (
            "phiRec",
            phi
        );
    }

    volVectorField&     HbyA    = *HbyAPtr;
    surfaceScalarField& rAUf    = *rAUfPtr;
    surfaceScalarField& phiFL   = *phiFLPtr;
    surfaceScalarField& phipB   = *phipBPtr;
    surfaceScalarField& phiHbyA = *phiHbyAPtr;
    surfaceScalarField& phiRec  = *phiRecPtr;

    // FIXME [Low]: Does not work with phiHbyA (only for twoFluids)
    // # include "scalePhi.H"

    // Non-orthogonal pressure corrector loop
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix* pEqnPtr = NULL;

        // Pressure corrector
        if (lorentzForceNonAlberto)
        {
            pEqnPtr = new fvScalarMatrix
            (
                fvm::laplacian(rAU, p)
             == fvc::div(phi)
            );
        }
        else
        {
            pEqnPtr = new fvScalarMatrix
            (
                fvm::laplacian(rAUf, p)
             == fvc::div(phiHbyA)
            );
        }

        fvScalarMatrix& pEqn = *pEqnPtr;

        // FIXME [Low]: Parallel mode: Reference setting not yet implemented (only for twoFluids)
        // # include "setReference.H"

        if
        (
            corr == nCorr-1
         && nonOrth == nNonOrthCorr
        )
        {
            pEqn.solve(mesh.solutionDict().solver("pFinal"));
        }
        else
        {
            pEqn.solve();
        }

        if (nonOrth == nNonOrthCorr)
        {
            // Flux corrector
            if (lorentzForceNonAlberto)
            {
                phi -= pEqn.flux();
            }
            else
            {
                phi = phiHbyA - pEqn.flux();

                phiRec = - pEqn.flux();
            }
        }

        // Delete p pointer
        delete pEqnPtr;
        pEqnPtr = NULL;
    }

#   include "continuityErrs.H"

    // Momentum corrector
    if (lorentzForceNonAlberto)
    {
        UcorrByrAU = -fvc::grad(p);

        // Explicitly relax UcorrByrAU for momentum corrector
        if (outerCorr != nOuterCorr - 1)
        {
            Info << "Relaxing UcorrByrAU with factor: " << mesh.solutionDict().relaxationFactor("p") << nl;
            UcorrByrAU.relax(mesh.solutionDict().relaxationFactor("p"));
        }

        U += rAU * UcorrByrAU;
    }
    else
    {
        phiRec += phiFL + phipB;

        UcorrByrAU = fvc::reconstruct( phiRec / rAUf );

        // Explicitly relax UcorrByrAU for momentum corrector
        if (outerCorr != nOuterCorr - 1)
        {
            Info << "Relaxing UcorrByrAU with factor: " << mesh.solutionDict().relaxationFactor("p") << nl;
            UcorrByrAU.relax(mesh.solutionDict().relaxationFactor("p"));
        }

        U = HbyA + rAU * UcorrByrAU;

        delete HbyAPtr;    HbyAPtr    = NULL;
        delete rAUfPtr;    rAUfPtr    = NULL;
        delete phiFLPtr;   phiFLPtr   = NULL;
        delete phipBPtr;   phipBPtr   = NULL;
        delete phiHbyAPtr; phiHbyAPtr = NULL;
        delete phiRecPtr;  phiRecPtr  = NULL;
    }

    // Correct boundary conditions
    U.correctBoundaryConditions();

//     Uphi = fvc::reconstruct(phi);
//     UmeshPhi = fvc::reconstruct(fvc::meshPhi(U));