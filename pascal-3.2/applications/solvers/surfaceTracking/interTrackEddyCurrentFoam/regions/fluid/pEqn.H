
// TODO [High]: Fix pressure relaxation for Alberto-method

    volScalarField rAU ("rAU", 1.0/UEqn.A());

    volVectorField*     HbyAPtr    = NULL;
    surfaceScalarField* rAUfPtr    = NULL;
    surfaceScalarField* phiFPtr    = NULL;
    surfaceScalarField* phiHbyAPtr = NULL;
    surfaceScalarField* phiRecPtr  = NULL;

    if (methodUp == "default")
    {
        U = rAU * UEqn.H();

        phi = (fvc::interpolate(U) & mesh.Sf());
    }
    else if (methodUp == "alberto")
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

        phiFPtr = new surfaceScalarField
        (
            "phiF",
            *rAUfPtr * ( fvc::interpolate(F) & mesh.Sf() )
        );

        phiHbyAPtr = new surfaceScalarField
        (
            "phiHbyA",
            ( fvc::interpolate(*HbyAPtr) & mesh.Sf() ) + *phiFPtr
        );

        phiRecPtr = new surfaceScalarField
        (
            "phiRec",
            phi
        );
    }

    volVectorField&     HbyA    = *HbyAPtr;
    surfaceScalarField& rAUf    = *rAUfPtr;
    surfaceScalarField& phiF    = *phiFPtr;
    surfaceScalarField& phiHbyA = *phiHbyAPtr;
    surfaceScalarField& phiRec  = *phiRecPtr;

// FIXME [Low]: Does not work with phiHbyA (only for twoFluids)
// # include "scalePhi.H"

    // Non-orthogonal pressure corrector loop
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix* pEqnPtr = NULL;

        // Pressure corrector
        if (methodUp == "default")
        {
            pEqnPtr = new fvScalarMatrix
            (
                fvm::laplacian(rAU, p)
             == fvc::div(phi)
            );
        }
        else if (methodUp == "alberto")
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
            if (methodUp == "default")
            {
                phi -= pEqn.flux();
            }
            else if (methodUp == "alberto")
            {
                phi = phiHbyA - pEqn.flux();

                phiRec = phiF - pEqn.flux();
            }
        }

        // Delete p pointer
        delete pEqnPtr;
        pEqnPtr = NULL;
    }

#   include "continuityErrs.H"

    // Momentum corrector
    if (methodUp == "default")
    {
//         UcorrByrAU = -fvc::grad(p);
//
//         // Explicitly relax UcorrByrAU for momentum corrector
//         if (outerCorr != nOuterCorr - 1)
//         {
//             Info << "Relaxing UcorrByrAU with factor: " << mesh.solutionDict().relaxationFactor("p") << nl;
//             UcorrByrAU.relax(mesh.solutionDict().relaxationFactor("p"));
//         }
//
//         U += rAU * UcorrByrAU;

        // Explicitly relax p for momentum corrector
        if (outerCorr != nOuterCorr - 1)
        {
            scalar alpha = mesh.solutionDict().relaxationFactor("p");

            if (alpha < 1.0)
            {
                Info << "Relaxing p with factor: "
                    << alpha << nl;

                p.relax(alpha);
            }
        }

        U -= rAU * fvc::grad(p);
    }
    else if (methodUp == "alberto")
    {
//         UcorrByrAU = fvc::reconstruct( phiRec / rAUf );
//
//         // Explicitly relax UcorrByrAU for momentum corrector
//         if (outerCorr != nOuterCorr - 1)
//         {
//             Info << "Relaxing UcorrByrAU with factor: " << mesh.solutionDict().relaxationFactor("p") << nl;
//             UcorrByrAU.relax(mesh.solutionDict().relaxationFactor("p"));
//         }
//
//         U = HbyA + rAU * UcorrByrAU;

        U = HbyA + rAU * fvc::reconstruct( phiRec / rAUf );

        delete HbyAPtr;    HbyAPtr    = NULL;
        delete rAUfPtr;    rAUfPtr    = NULL;
        delete phiFPtr;    phiFPtr    = NULL;
        delete phiHbyAPtr; phiHbyAPtr = NULL;
        delete phiRecPtr;  phiRecPtr  = NULL;
    }

    U.correctBoundaryConditions();

    if (debug)
    {
        volScalarField& debug_divU = *debug_divUPtr;
        volScalarField& debug_divPhi = *debug_divPhiPtr;
        volVectorField& debug_Uphi = *debug_UphiPtr;
        volVectorField& debug_UmeshPhi = *debug_UmeshPhiPtr;

        scalar meshH = pow(max(mesh.V()),1.0/3.0).value();

        debug_divU = fvc::div(U) * meshH;
        debug_divPhi = fvc::div(phi) * meshH;

        debug_Uphi = fvc::reconstruct(phi);
        debug_UmeshPhi = fvc::reconstruct(fvc::meshPhi(U));
    }
