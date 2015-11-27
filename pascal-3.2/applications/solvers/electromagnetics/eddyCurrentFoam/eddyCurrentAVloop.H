
# define updateSigmaA(regionI)                                                \
                                                                              \
    sigma_.mapExtrapolate(regionI);                                           \
    ARe_.mapInterpolate(regionI);                                             \
    AIm_.mapInterpolate(regionI);

# define updateRmapGradient(regionI,vsf)                                      \
                                                                              \
    vsf##Grad_[regionI] = fvc::grad(vsf);                                     \
    vsf##Grad_.rmap(regionI);

# define updateGradientV(regionI)                                             \
                                                                              \
    updateRmapGradient(regionI,VRe);                                          \
    updateRmapGradient(regionI,VIm);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // AV iterations
    while (control.loop())
    {
        // Init gradient of V in first iteration
        if (control.firstIter())
        {
            updateGradientV(control.conductor());
        }

        // Solve for A in base region
        {
            setRegionScope(control.base());

            // Relax gradient of V
            control.relax(VReGrad);
            control.relax(VImGrad);

#           include "eddyCurrentAEqn.H"
        }

        // Update A and sigma
        updateSigmaA(control.conductor());

        // Solve for V in conductor region
        {
            setRegionScope(control.conductor());

#           include "eddyCurrentVEqn.H"
        }

        // Update gradient of V
        updateGradientV(control.conductor());

        control.subWrite();
    }
