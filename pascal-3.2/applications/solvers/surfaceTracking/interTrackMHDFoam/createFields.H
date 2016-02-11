
    Info << nl;

    Info << "Reading field p" << nl << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info << "Reading field U" << nl << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    #   include "createPhi.H"


//     Info << "Creating field Uphi" << nl << endl;
//     volVectorField Uphi
//     (
//         IOobject
//         (
//             "Uphi",
//             runTime.timeName(),
//             mesh,
//             IOobject::NO_READ,
//             IOobject::AUTO_WRITE
//         ),
//         fvc::reconstruct(phi)
//     );
//     Uphi.write();
//
//
//     Info << "Creating field UmeshPhi" << nl << endl;
//     volVectorField UmeshPhi
//     (
//         IOobject
//         (
//             "UmeshPhi",
//             runTime.timeName(),
//             mesh,
//             IOobject::NO_READ,
//             IOobject::AUTO_WRITE
//         ),
//         fvc::reconstruct(fvc::meshPhi(U))
//     );
//     UmeshPhi.write();


    Info << "Reading field FL" << nl << endl;
    volVectorField FL
    (
        IOobject
        (
            "FL",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info << "Reading field pB" << nl << endl;
    volScalarField pB
    (
        IOobject
        (
            "pB",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    volScalarField* pBfacPtr = NULL;
    {
#       include "readControls.H"

        pBfacPtr = new volScalarField
        (
            "pBfac",
            pB * lorentzForceRotationalFactor
        );
    }
    volScalarField& pBfac = *pBfacPtr;


    Info << "Reading field fluidIndicator" << nl << endl;
    volScalarField fluidIndicator
    (
        IOobject
        (
            "fluidIndicator",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info << "Reading g" << nl << endl;
    uniformDimensionedVectorField g
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    Info << "Reading transportProperties" << endl;
    twoPhaseMixture twoPhaseProperties(U, phi, "fluidIndicator");


    // Construct incompressible turbulence model
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, twoPhaseProperties)
    );
    Info << nl;


    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        twoPhaseProperties.rho()
    );
    surfaceScalarField rhof
    (
        IOobject
        (
            "rhof",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rho)
    );


    volScalarField mu
    (
        IOobject
        (
            "mu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        twoPhaseProperties.mu()
    );

    {
#       include "readControls.H"
#       include "updateLorentzForcePressureCorr.H"
    }


    Info << "Creating tracked surface" << nl << endl;
    trackedSurface interface(mesh, rho, U, p, pBfac, phi, g, twoPhaseProperties, turbulence);
    Info << nl;


    // Read interpolators if present
    if (interface.twoFluids())
    {
        interface.interpolatorAB();
    }


    // FIXME [Low]: Parallel mode: Reference setting not yet implemented (only for twoFluids)
    // #   include "setRefCell.H"



    Info << "Creating field field UcorrByrAU" << nl << endl;
    volVectorField UcorrByrAU
    (
        IOobject
        (
            "UcorrByrAU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::grad(p)
    );


    Info << "Creating/calculating field rAU" << endl;
    volScalarField* rAUPtr = NULL;
    {
#       include "readControls.H"

        int outerCorr = nOuterCorr-1;

#       include "UEqn.H"

        rAUPtr = new volScalarField
        (
            IOobject
            (
                "rAU",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            1.0/UEqn.A()
        );

        Info << nl;
    }
    volScalarField& rAU = *rAUPtr;