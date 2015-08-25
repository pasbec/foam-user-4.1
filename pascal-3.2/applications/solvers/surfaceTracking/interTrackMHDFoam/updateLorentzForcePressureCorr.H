
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

            // Normal vectors
            vectorField nfp = mesh.boundary()[patchI].nf();

            // Deltas
            vectorField deltap = mesh.boundary()[patchI].delta();

            // DeltaCoeffs
            scalarField deltaCoeffsp = mesh.boundary()[patchI].deltaCoeffs();

            // Lorentz force
            vectorField FLp = FL.boundaryField()[patchI];

            // Magnetic pressure gradient
            scalarField pBsnGradp (FLp.size(), 0.0);

            if (lorentzForceRotationalFactor > 0.0)
            {
                if
                (
                    Foam::word
                    (
                        mesh.schemesDict().gradScheme
                        (
                            "grad(pB)"
                        )
                    ) != "extendedLeastSquares"
                )
                {
                    Warning << "The gradScheme for pB "
                            << "should be \"extendedLeastSquares 0\" for the boundary "
                            << "non-orthogonal correction to be right" << endl;
                }

                // Gradient field of pB with non-orthogonal correction at boundaries
                volVectorField pBGrad = fvc::grad(pB);

                // Correction vector
                vectorField kp = deltap - nfp*(nfp&deltap);

                // Magnetic pressure sn-gradient with non-orthogonal correction at boundaries
                pBsnGradp = (
                    pB.boundaryField()[patchI]
                  - pB.boundaryField()[patchI].patchInternalField()
                  - (kp & pBGrad.boundaryField()[patchI].patchInternalField())
                ) * deltaCoeffsp;

//                 // Magnetic pressure gradient without correction
//                 scalarField pBsnGradp = pB.boundaryField()[patchI].snGrad();
            }

            // Density
            scalarField rhop = rho.boundaryField()[patchI];

            // Set correct pressure gradient
            pp.gradient() = rhop * (FLp & nfp);
            pp.gradient() += pBsnGradp * lorentzForceRotationalFactor;
            pp.gradient() *= lorentzForceVolumeFactor;
        }
    }