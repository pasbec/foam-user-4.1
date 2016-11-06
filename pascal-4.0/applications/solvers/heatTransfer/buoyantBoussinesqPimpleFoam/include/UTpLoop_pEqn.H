
    {
        volScalarField rAU("rAU", 1.0/UEqn.A());
        surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

        U = rAU*UEqn.H();

        surfaceScalarField phiU
        (
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rAU, U, phi)
        );

        phi = phiU + rAUf*fvc::interpolate(rhok)*(g & mesh.Sf());

        // Manually correct the boundary conditions for the pressure at walls
        // if the additional volume forcing is activated. For this, it is
        // utterly important to use fixedGradient patchFields instead of
        // the zeroGradient ones at all walls. bouyantPressure may not be used, too
// TODO: Do not use fixedGradient
        if
        (
            storage.item_F().enabled()
        )
        {
            forAll(p.boundaryField(), patchI)
            {
                if (isA<wallFvPatch>(mesh.boundary()[patchI]))
                {
                    if
                    (
                        p.boundaryField()[patchI].type()
                     == fixedFluxPressureFvPatchScalarField::typeName
                    )
                    {
                        if (Control::debug)
                        {
                            Info<< Control::typeName << " | pEqn.H : "
                                << "Pressure at patch "
                                << mesh.boundary()[patchI].name()
                                << "is being corrected by means of a "
                                << "'fixedFluxPressure' boundary condition. "
                                << "This is OK."
                                << endl;
                        }
                    }
                    else if
                    (
                        p.boundaryField()[patchI].type()
                     == fixedGradientFvPatchScalarField::typeName
                    )
                    {
                        fixedGradientFvPatchScalarField& pp =
                            refCast<fixedGradientFvPatchScalarField>
                            (
                                p.boundaryField()[patchI]
                            );

                        // Patch data
                        const vectorField& Sfp =
                            mesh.Sf().boundaryField()[patchI];
                        const scalarField& magSfp =
                            mesh.magSf().boundaryField()[patchI];


                        if (!settings.snGradpFromFlux)
                        {
                            if (Control::debug)
                            {
                                Info<< Control::typeName << " | pEqn.H : "
                                    << "Correcting pressure at patch "
                                    << mesh.boundary()[patchI].name() << "."
                                    << endl;
                            }

                            // Density
                            const scalarField& rhokp =
                                rhok.boundaryField()[patchI];

                            // Direct way by means of F
                            const vectorField& Fp =
                                storage.F().boundaryField()[patchI];

                            pp.gradient() = ((rhokp*g.value() + Fp) & Sfp)/magSfp;
                        }
                        else
                        {
                            if (Control::debug)
                            {
                                Info<< Control::typeName << " | pEqn.H : "
                                    << "Correcting pressure at patch "
                                    << mesh.boundary()[patchI].name()
                                    << " from flux-difference."
                                    << endl;
                            }

                            // Indirect way by means of fixed flux-pressure
                            const scalarField& phip =
                                phi.boundaryField()[patchI];
                            const vectorField& Up =
                                U.boundaryField()[patchI];
                            const scalarField& rAUfp =
                                rAUf.boundaryField()[patchI];

                            pp.gradient() = ((phip - (Up & Sfp))/rAUfp)/magSfp;
                        }

                        pp.updateCoeffs();
                    }
                    else
                    {
                        WarningIn(Control::typeName + " | pEqn.H")
                            << "The pressure boundary condition for "
                            << "wall-patch named '"
                            << mesh.boundary()[patchI].name() << "' is '"
                            << p.boundaryField()[patchI].type() << "',"
                            << endl << "    "
                            << "It should be of type 'fixedGradient' "
                            << "or 'fixedFluxPressure' "
                            << "to correctly represent the gradient dp/dn."
                            << endl;
                    }
                }
            }

            p.correctBoundaryConditions();
        }

        while (control.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian(rAUf, p) == fvc::div(phi)
            );

// TODO: Read pRefCell and pRefValue with control
// #           include "setReference.H"
//             pEqn.setReference(pRefCell, pRefValue);

            pEqn.solve(mesh.solutionDict().solver(p.select(control.finalInnerIter())));

            if (control.finalNonOrthogonalIter())
            {
                phi -= pEqn.flux();
            }
        }

        U += rAU*fvc::reconstruct((phi - phiU)/rAUf);
        U.correctBoundaryConditions();

// TODO: Continuity errors
// #       include "continuityErrs.H"
    }
