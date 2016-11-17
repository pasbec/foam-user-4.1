
    {
        alphat = turbulence.nut()/Prt;
        alphat.correctBoundaryConditions();

        volScalarField kappaEff("kappaEff", turbulence.nu()/Pr + alphat);

        fvScalarMatrix TEqn
        (
            fvm::ddt(T)
          + fvm::div(phi, T)
          - fvm::laplacian(kappaEff, T)
//          ==
//             radiation.ST(rhoCpRef, T)
        );

        TEqn.relax();

        TEqn.solve();

    //     radiation.correct();
    }
