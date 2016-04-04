
    volScalarField rAU(1.0/UEqn().A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn().H();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh.Sf())
// TODO: Port ddtCorr from OpenFOAM 3.0.x
//     + rAUf*fvc::ddtCorr(U, phi)
    );

    autoPtr<surfaceScalarField> phiF;

    if
    (
        storage.is_F()
    && !settings.UpDirectForce
    )
    {
        if (globalSettings.debug)
        {
            Info << "interTrackManager(pEqn.H) : "
                << "Adding volume force flux to phiHbyA."
                    << endl;
        }

        phiF.set
        (
            new surfaceScalarField
            (
                "phiF",
                rAUf*(fvc::interpolate(storage.F()) & mesh.Sf())
            )
        );

        phiHbyA += phiF();
    }

// TODO: Two fluids
// #   include "scalePhi.H"
// if (p.needReference())
// {
//     fvc::makeRelative(phiHbyA, U);
//     adjustPhi(phiHbyA, U, p);
//     fvc::makeAbsolute(phiHbyA, U);
// }

    tmp<volScalarField> rAtU(rAU);

// TODO: Learn and get to know if this is ok with F
//       but I have the feeling, that it is NOT!
    if (control.consistent())
    {
        if (globalSettings.debug)
        {
            Info << "interTrackManager(pEqn.H) : "
                << "Applying consistent PIMPLE correction."
                    << endl;
        }

        rAtU = 1.0/max(1.0/rAU - UEqn().H1(), 0.1/rAU);
        phiHbyA +=
            fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
        HbyA -= (rAU - rAtU())*fvc::grad(p);
    }

    rAUf = fvc::interpolate(rAtU());

    if (control.nCorrPISO() <= 1)
    {
        UEqn.clear();
    }

    // Manually correct the boundary conditions for the pressure at walls
    // if the addititional volume forcing is activated. For this, it is
    // utterly important to use fixedGradient patchFields instead of
    // the zeroGradient ones at all walls
// TODO: What about other pressure boundaries? It is obviously better
//       to conduct a more general aproach. Therefor it makes probably
//       sense to port back the fixedFluxPressureFvPatchFields from
//       OpenFOAM 3.0.x, does it!?
    if (storage.is_F())
    {
        forAll(p.boundaryField(), patchI)
        {
            if
            (
                p.boundaryField()[patchI].type()
            == fixedGradientFvPatchField<scalar>::typeName
            )
            {
                fixedGradientFvPatchField<scalar>& pp =
                    refCast<fixedGradientFvPatchField<scalar> >
                    (
                        p.boundaryField()[patchI]
                    );

                // Patch data
                const vectorField& Sfp =
                    mesh.Sf().boundaryField()[patchI];
                const scalarField& mSfp =
                    mesh.magSf().boundaryField()[patchI];


                if (settings.pFdirectCorrection)
                {
                    if (globalSettings.debug)
                    {
                        Info << "interTrackManager(pEqn.H) : "
                            << "Correcting pressure at patch "
                                << mesh.boundary()[patchI].name()
                                << " directly, by means of volume force."
                                << endl;
                    }

                    // Direct way by means of F
                    const vectorField& Fp =
                        storage.F().boundaryField()[patchI];

                    pp.gradient() = (Fp & Sfp)/mSfp;
                }
                else
                {
                    if (globalSettings.debug)
                    {
                        Info << "interTrackManager(pEqn.H) : "
                            << "Correcting pressure at patch "
                                << mesh.boundary()[patchI].name()
                                << " indirectly, by means of fixed flux-pressure."
                                << endl;
                    }

                    // Indirect way by means of fixed flux-pressure
                    const scalarField& phiHbyAp =
                        phiHbyA.boundaryField()[patchI];
                    const vectorField& Up =
                        U.boundaryField()[patchI];
                    const scalarField& rAUfp =
                        rAUf.boundaryField()[patchI];

                    pp.gradient() = rAUfp*(phiHbyAp - (Up & Sfp))/mSfp;
                }
            }
        }
    }

    // Non-orthogonal pressure corrector loop
    while (control.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAUf, p) == fvc::div(phiHbyA)
        );

// TODO: Two fluids
// TODO: Read pRefCell and pRefValue with control
// #include "setReference.H"
//         pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve
        (
            mesh.solutionDict().solverDict
            (
                p.select(control.finalInnerIter())
            )
        );

        // Manually store residuls after first
        // non-orthogonal iteration
        if (control.firstInnerIter())
        {
            control.storeResiduals(p.name());
        }

        // Stop non-orthogonal correction after the next (final)
        // iteration, if the last number of solution iterations
        // was "0" for the named field
        control.skipZeroNonOrtho(p.name());

        if (control.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p.relax();

            U = HbyA;

            if
            (
                storage.is_F()
            && !settings.UpDirectForce
            )
            {
                // Correct the momentum source with the pressure gradient
                // flux calculated from the relaxed pressure
                U += rAtU*fvc::reconstruct((phiF() - pEqn.flux())/rAUf);

                phiF.clear();
            }
            else
            {
                // Correct the momentum source with the pressure gradient
                // calculated from the relaxed pressure, directly
                U -= rAtU*fvc::grad(p);
            }

            U.correctBoundaryConditions();
        }
    }

#   include "continuityErrs.H"

// TODO: Port ddtCorr from OpenFOAM 3.0.x
//     {
//         Uf = fvc::interpolate(U);
//         surfaceVectorField n(mesh.Sf()/mesh.magSf());
//         Uf += n*(phi/mesh.magSf() - (n & Uf));
//     }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

