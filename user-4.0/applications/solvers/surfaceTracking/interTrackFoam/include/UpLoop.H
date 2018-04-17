

    {
        if (Control::debug)
        {
            Info<< Control::typeName << " | UpLoop.H : "
                << "Commencing PIMPLE U-p loop."
                << endl;
        }

        // --- PIMPLE corrector loop
        while (control.loop())
        {
            uniformDimensionedVectorField& g = storage.g();
            volScalarField& p = storage.p();
            volVectorField& U = storage.U();
            surfaceScalarField& phi = storage.phi();
// TODO: Make trackedSurface use real rho and mu fields!
            volScalarField& rho = storage.rho();
            twoPhaseMixture& transport = storage.transport();
            incompressible::turbulenceModel& turbulence = storage.turbulence();
            trackedSurface& interface = storage.interface();

            if (settings.relToUinf)
            {
                // Update free-stream force
// TODO: Check force calculation!
                storage.fInf().value() =
                    interface.totalPressureForce();
//                 storage.fInf().value() =
//                   - interface.totalPressureForce()
//                   + interface.totalViscousForce();
//                 storage.fInf().value() =
//                     interface.totalPressureForce()
//                   + interface.normalViscousForce()
//                   + interface.tangentialViscousForce();

                dimensionedVector Finf
                (
                    "Finf",
                    dimForce,
                    vector::zero
                );

                forAll(mesh.boundary(), patchI)
                {
                    if
                    (
                        mesh.boundary()[patchI].type()
                     != wedgeFvPatch::typeName
                    )
                    {
                        const vectorField& Sf = mesh.Sf().boundaryField()[patchI];

                        const scalarField& P = p.boundaryField()[patchI];

                        symmTensorField tau =
                            (
                                rho*turbulence.nu()
                            * dev(twoSymm(fvc::grad(U)))
                            )().boundaryField()[patchI];
//                           - rho.boundaryField()[patchI]
//                           * turbulence.devReff()().boundaryField()[patchI];

                        Finf.value() += gSum(Sf & (P*I - tau));
                    }
                }

                dimensionedVector FinfVol
                (
                    "Finf",
                    fvc::domainIntegrate
                    (
                        fvc::div
                        (
                            p*I
                          - rho*turbulence.nu()*dev(twoSymm(fvc::grad(U)))
                        )
                    )
                );

                // Constraint free-stream force
                forAll(mesh.boundary(), patchI)
                {
                    if
                    (
                        mesh.boundary()[patchI].type()
                     == wedgeFvPatch::typeName
                    )
                    {
                        const wedgePolyPatch& wp =
                            refCast<const wedgePolyPatch>
                            (
                                mesh.boundary()[patchI].patch()
                            );

                        vector ez = wp.axis();
                        vector ephi = wp.centreNormal();
                        vector er = ez ^ ephi;

//                         storage.fInf().value() -=
//                             sqr(er) & storage.fInf().value();

                        storage.fInf().value() -=
                            sqr(vector(1,0,0)) & storage.fInf().value();
                        storage.fInf().value() -=
                            sqr(vector(0,0,1)) & storage.fInf().value();

                        Finf.value() -=
                            sqr(vector(1,0,0)) & Finf.value();
                        Finf.value() -=
                            sqr(vector(0,0,1)) & Finf.value();

                        FinfVol.value() -=
                            sqr(vector(1,0,0)) & FinfVol.value();
                        FinfVol.value() -=
                            sqr(vector(0,0,1)) & FinfVol.value();
                    }
                }

                // Store old free-stream velocity
                dimensionedVector Uinf(storage.Uinf());
                if (control.firstIter())
                {
                    storage.Uinf_0() = storage.Uinf();
                }

                // Calculate acceleration from du = F/m * dt
                dimensionedVector aInf
                (
                    "aInf",
                    storage.fInf() / storage.mInf()

                );

//                 // Calculate acceleration from du = F/m * dt
//                 dimensionedVector aInf
//                 (
//                     "aInf",
//                     FinfVol / fvc::domainIntegrate(rho)
//
//                 );

                // Calculate change of Uinf
                dimensionedVector dUinf
                (
                    "dUinf",
                   - aInf * runTime.deltaT()
                );

                // Update free-stream velocity with relaxation
                storage.Uinf() =
                    storage.Uinf_0() + storage.dUinfRelax() * dUinf;

                Info<< "Free-stream velocity adjustment:"
                    << " fInf = " << storage.fInf().value()
                    << " , Uinf = " << storage.Uinf().value()
                    << " , Uinf_0 = " << storage.Uinf_0().value()
                    << " , deltaT = " << runTime.deltaT().value()
                    << " , aInf = " << aInf.value()
                    << " , dUinf = " << dUinf.value()
                    << " , max(mag(U)) = " << max(mag(U)).value()
                    << " , max(mag(U0)) = " << max(mag(storage.U0())).value()
                    << " , Finf = " << Finf.value()
                    << " , FinfVol = " << FinfVol.value()
                    << endl;

                // Adjust boundary conditions for velocity
                {
                    forAll(U.boundaryField(), patchI)
                    {
                        if
                        (
                            U.boundaryField()[patchI].type()
                        == fixedValueFvPatchVectorField::typeName
                        && patchI != interface.aPatchID()
                        )
                        {
                            U.boundaryField()[patchI] ==
                                storage.Uinf().value();
                        }
                    }

                    U.correctBoundaryConditions();
                }

//                 // Remove old counterflow velocity from U
//                 {
//                     // WARNING: FORCE assign boundary values!
//                     U == U - storage.U0();
//                     U.correctBoundaryConditions();
//
//                     phi -= storage.phi0();
//                 }
//
//                 // Adjust boundary conditions for counterflow velocity
//                 {
//                     storage.U0() *= mag(storage.Uinf().value())/(mag(Uinf.value()) + VSMALL);
//                     storage.phi0() *= mag(storage.Uinf().value())/(mag(Uinf.value()) + VSMALL);
//
//                     forAll(storage.U0().boundaryField(), patchI)
//                     {
//                         if
//                         (
//                             storage.U0().boundaryField()[patchI].type()
//                         == fixedValueFvPatchVectorField::typeName
//                         && patchI != interface.aPatchID()
//                         )
//                         {
//                             storage.U0().boundaryField()[patchI] ==
//                                 storage.Uinf().value();
//                         }
//                     }
//
//                     storage.U0().correctBoundaryConditions();
//                 }
//
//                 // Update counterflow velocity
//                 {
//                     volScalarField& p0 = storage.p0();
//                     volVectorField& U0 = storage.U0();
//                     surfaceScalarField& phi0 = storage.phi0();
//
//                     volScalarField mu
//                     (
//                         "mu",
//                         rho*transport.nu()
//                     );
//
//                     surfaceScalarField muf
//                     (
//                         "muf",
//                         fvc::interpolate(mu)
//                     );
//
//                     surfaceScalarField rhoPhi0
//                     (
//                         "rho*phi0",
//                         fvc::interpolate(rho)*phi0
//                     );
//
//                     fvVectorMatrix U0Eqn
//                     (
//                         fvm::ddt(rho, U0)
//                       + fvm::div(rhoPhi0, U0)
//                       - fvm::laplacian(muf, U0)
//                       - (fvc::grad(U0) & fvc::grad(muf))
//                     );
//
//                     U0Eqn.relax();
//
//                     solve(U0Eqn == -fvc::grad(p0));
//
//                     while (control.correct())
//                     {
//                         volScalarField rAU("rAU", 1.0/U0Eqn.A());
//                         surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
//
//                         U0 = rAU*U0Eqn.H();
//
//                         phi0 = fvc::interpolate(U0) & mesh.Sf();
//
//                         // Non-orthogonal pressure corrector loop
//                         while (control.correctNonOrthogonal())
//                         {
//                             fvScalarMatrix p0Eqn
//                             (
//                                 fvm::laplacian(rAUf,p0) == fvc::div(phi0)
//                             );
//
//                             p0Eqn.setReference
//                             (
//                                 control.pRefCell(),
//                                 control.pRefValue()
//                             );
//
//                             p0Eqn.solve();
//
//                             if (control.finalNonOrthogonalIter())
//                             {
//                                 phi0 -= p0Eqn.flux();
//                             }
//                         }
//
//                         p0.relax();
//
//                         U0 -= rAU*fvc::grad(p0);
//                         U0.correctBoundaryConditions();
//                     }
//                 }
//
//                 // Add new counterflow velocity to U
//                 {
//                     // WARNING: FORCE assign boundary values!
//                     U == U + storage.U0();
//                     U.correctBoundaryConditions();
//
//                     phi += storage.phi0();
//                 }

                // WARNING: WE DO NOT TOUCH boundary values!
                storage.c() += storage.gradcInf() & mesh.C();
            }

            interface.updateBoundaryConditions(); // abs phi

            if (settings.relToUinf)
            {
                // WARNING: WE DO NOT TOUCH boundary values!
                storage.c() -= storage.gradcInf() & mesh.C();
            }

            surfaceScalarField gxf("gxf", g & mesh.Cf());

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phi, U);

#           include "UpLoop_UEqn.H"

            // --- Pressure corrector PISO loop
            while (control.correct())
            {
#               include "UpLoop_pEqn.H"
            }

            if (control.turbCorr())
            {
                turbulence.correct();
            }

#           include "UpLoop_cEqn.H"
#           include "UpLoop_rhoUpdate.H"

            // Make the fluxes absolute
            fvc::makeAbsolute(phi, U);

            interface.updatePoints(); // abs phi

// TODO: Continuity errors
// #           include "trackedSurfaceContinuityErrs.H"
        }
    }
