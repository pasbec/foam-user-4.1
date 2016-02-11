
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
            tmp<vectorField> tnfp
            (
                new vectorField
                (
                    mesh.boundary()[patchI].nf()
                )
            );

            vectorField& nfp = tnfp();

            // Volume force
            const vectorField& Fp = F.boundaryField()[patchI];

            // Set correct pressure gradient
            pp.gradient() = (Fp & nfp);

            tnfp.clear();
        }
    }
