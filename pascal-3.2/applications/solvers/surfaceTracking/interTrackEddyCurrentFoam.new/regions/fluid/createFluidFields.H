
    // Gravity
    Info << "Create field g for"
        << " region " << mesh_[control.fluidA()].name()
        << nl << endl;
    uniformDimensionedVectorField g
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh_[control.fluidA()],
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    // Pressure
    Info << "Create field p for"
        << " region " << mesh_[control.fluidA()].name()
        << nl << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh_[control.fluidA()],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_[control.fluidA()]
    );


    // Velocity
    Info << "Create field U for"
        << " region " << mesh_[control.fluidA()].name()
        << nl << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh_[control.fluidA()],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_[control.fluidA()]
    );


    volVector4Field* UpPtr = NULL;
    {
#       include "readFluidControls.H"

        if (methodUp == "coupled")
        {
            Info << "Create field Up for"
                << " region " << mesh_[control.fluidA()].name()
                << nl << endl;
            UpPtr = new volVector4Field
            (
                IOobject
                (
                    "Up",
                    runTime.timeName(),
                    mesh_[control.fluidA()],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_[control.fluidA()],
                dimensionedVector4("zero", dimless, vector4::zero)
            );
        }
    }


    // Face volume flux based on U
    Info << "Create field phi for"
        << " region " << mesh_[control.fluidA()].name()
        << nl << endl;
    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh_[control.fluidA()],
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U) & mesh_[control.fluidA()].Sf()
    );


    // Phase indicator
    Info << "Create field fluidIndicator for"
        << " region " << mesh_[control.fluidA()].name()
        << nl << endl;
    volScalarField fluidIndicator
    (
        IOobject
        (
            "fluidIndicator",
            runTime.timeName(),
            mesh_[control.fluidA()],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_[control.fluidA()]
    );


    // Transport model
    Info << "Create transport model for"
        << " region " << mesh_[control.fluidA()].name()
        << endl;
    twoPhaseMixture twoPhaseProperties(U, phi, "fluidIndicator");
    Info << endl;


    // Density
    Info << "Create field rho for"
        << " region " << mesh_[control.fluidA()].name()
        << nl << endl;
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh_[control.fluidA()],
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        twoPhaseProperties.rho()
    );


    // Dynamic viscosity
    Info << "Create field mu for"
        << " region " << mesh_[control.fluidA()].name()
        << nl << endl;
    volScalarField mu
    (
        IOobject
        (
            "mu",
            runTime.timeName(),
            mesh_[control.fluidA()],
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        twoPhaseProperties.mu()
    );


    // Turbulence model
    Info << "Create turbulence model for"
        << " region " << mesh_[control.fluidA()].name()
        << endl;
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New
        (
            U,
            phi,
            twoPhaseProperties
        )
    );
    Info << nl;


//     // Tracked surface
//     Info << "Create tracked surface for"
//         << " region " << mesh_[control.fluidA()].name()
//         << endl;
//     trackedSurface interface
//     (
//         mesh_[control.fluidA()],
//         rho,
//         U,
//         p,
//         pB_[control.fluidA()],
//         phi,
//         g,
//         twoPhaseProperties,
//         turbulence
//     );
//     Info << nl;

    // Tracked surface
    Info << "Create tracked surface for"
        << " region " << mesh_[control.fluidA()].name()
        << endl;
    trackedSurface interface
    (
        mesh_[control.fluidA()],
        rho,
        U,
        p,
        phi,
        g,
        twoPhaseProperties
    );
    Info << nl;


    // Volume force
    Info << "Create field F for"
        << " region " << mesh_[control.fluidA()].name()
        << nl << endl;
    volVectorField F
    (
        IOobject
        (
            "F",
            runTime.timeName(),
            mesh_[control.fluidA()],
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_[control.fluidA()],
        dimensioned<vector>
        (
            "F",
            dimMass/pow(dimLength,2)/pow(dimTime,2),
            pTraits<vector>::zero
        ),
        calculatedFvPatchField<vector>::typeName
    );
    F.write();


    volVectorField* debug_FgradPtr = NULL;
    volScalarField* debug_divUPtr = NULL;
    volScalarField* debug_divPhiPtr = NULL;
    volVectorField* debug_UphiPtr = NULL;
    volVectorField* debug_UmeshPhiPtr = NULL;

    if (debug)
    {
        Info << "Create field Fgrad for"
            << " region " << mesh_[control.fluidA()].name()
            << nl << endl;
        debug_FgradPtr = new volVectorField
        (
            IOobject
            (
                "Fgrad",
                runTime.timeName(),
                mesh_[control.fluidA()],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                true
            ),
            -fvc::grad(pB_[control.fluidA()])
        );
        debug_FgradPtr->write();


        Info << "Create field divU for"
            << " region " << mesh_[control.fluidA()].name()
            << nl << endl;
        debug_divUPtr = new volScalarField
        (
            IOobject
            (
                "divU",
                runTime.timeName(),
                mesh_[control.fluidA()],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                true
            ),
            fvc::div(U)
        );
        debug_divUPtr->write();


        Info << "Create field divPhi for"
            << " region " << mesh_[control.fluidA()].name()
            << nl << endl;
        debug_divPhiPtr = new volScalarField
        (
            IOobject
            (
                "divPhi",
                runTime.timeName(),
                mesh_[control.fluidA()],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                true
            ),
            fvc::div(phi)
        );
        debug_divPhiPtr->write();


        Info << "Create field Uphi for"
            << " region " << mesh_[control.fluidA()].name()
            << nl << endl;
        debug_UphiPtr = new volVectorField
        (
            IOobject
            (
                "Uphi",
                runTime.timeName(),
                mesh_[control.fluidA()],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::reconstruct(phi)
        );
        debug_UphiPtr->write();


        Info << "Create field UmeshPhi for"
            << " region " << mesh_[control.fluidA()].name()
            << nl << endl;
        debug_UmeshPhiPtr = new volVectorField
        (
            IOobject
            (
                "UmeshPhi",
                runTime.timeName(),
                mesh_[control.fluidA()],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::reconstruct(fvc::meshPhi(U))
        );
        debug_UmeshPhiPtr->write();
    }
