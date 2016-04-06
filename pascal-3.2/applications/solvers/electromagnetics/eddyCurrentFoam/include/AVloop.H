
# define updateSigmaA(regionI)                                                \
                                                                              \
    sigma_.mapExtrapolate(regionI);                                           \
    ARe_.mapInterpolate(regionI);                                             \
    AIm_.mapInterpolate(regionI);

# define updateRmapGradient(regionI,vsf)                                      \
                                                                              \
    vsf##Grad_[regionI] == fvc::grad(vsf);                                    \
    vsf##Grad_.rmap(regionI);

# define updateGradientV(regionI)                                             \
                                                                              \
    updateRmapGradient(regionI,VRe);                                          \
    updateRmapGradient(regionI,VIm);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Init gradient of V in first iteration
    if (control.meshIs3D() )
    {
        updateGradientV(control.conductor());
    }

    // AV iterations
    while (control.loop())
    {
        // Solve for A in base region
        {
            setRegionScope(control.base());

#           include "AEqn.H"
        }

        // Only if mesh is 3D
        if (control.meshIs3D())
        {
            // Update A and sigma
            updateSigmaA(control.conductor());

            // Solve for V in conductor region
            {
                setRegionScope(control.conductor());

                // Relax A
                control.relax(ARe);
                control.relax(AIm);

#               include "VEqn.H"
            }

            // Update gradient of V
            updateGradientV(control.conductor());
        }
    }

    Info << nl;
