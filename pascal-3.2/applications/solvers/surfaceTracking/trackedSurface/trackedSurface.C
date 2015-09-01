/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "trackedSurface.H"

#include "volFields.H"
#include "transformField.H"

#include "emptyFaPatch.H"
#include "wedgeFaPatch.H"
#include "wallFvPatch.H"

#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"

#include "tetFemMatrices.H"
#include "tetPointFields.H"
#include "faceTetPolyPatch.H"
#include "tetPolyPatchInterpolation.H"
#include "fixedValueTetPolyPatchFields.H"
#include "fixedValuePointPatchFields.H"
#include "twoDPointCorrector.H"

#include "slipFvPatchFields.H"
#include "symmetryFvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "zeroGradientCorrectedFvPatchFields.H"
#include "fixedGradientCorrectedFvPatchFields.H"
#include "fixedValueCorrectedFvPatchFields.H"

#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(trackedSurface, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void trackedSurface::clearOut()
{
    deleteDemandDrivenData(interpolatorABPtr_);
    deleteDemandDrivenData(interpolatorBAPtr_);
    deleteDemandDrivenData(controlPointsPtr_);
    deleteDemandDrivenData(motionPointsMaskPtr_);
    deleteDemandDrivenData(pointsDisplacementDirPtr_);
    deleteDemandDrivenData(facesDisplacementDirPtr_);
    deleteDemandDrivenData(totalDisplacementPtr_);
    deleteDemandDrivenData(aMeshPtr_);
    deleteDemandDrivenData(UsPtr_);
    deleteDemandDrivenData(phisPtr_);
    deleteDemandDrivenData(surfactConcPtr_);
    deleteDemandDrivenData(surfaceTensionPtr_);
    deleteDemandDrivenData(surfactantPtr_);
    deleteDemandDrivenData(fluidIndicatorPtr_);
    deleteDemandDrivenData(muEffFluidAvalPtr_);
    deleteDemandDrivenData(muEffFluidBvalPtr_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

trackedSurface::trackedSurface
(
    dynamicFvMesh& m,
    const volScalarField& rho,
    volVectorField& Ub,
    volScalarField& Pb,
    const volScalarField& PbExt,
    const surfaceScalarField& sfPhi,
    const uniformDimensionedVectorField g,
    const twoPhaseMixture& transportModel,
    const autoPtr<incompressible::turbulenceModel>& turbulenceModelPtr
)
:
    IOdictionary
    (
        IOobject
        (
            "trackedSurfaceProperties",
            Ub.mesh().time().constant(),
            Ub.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(m),
    rho_(rho),
    U_(Ub),
    p_(Pb),
    pExt_(PbExt),
    phi_(sfPhi),
    g_(g),
    transport_(transportModel),
    turbulence_(turbulenceModelPtr),
    curTimeIndex_(Ub.mesh().time().timeIndex()),
    twoFluids_
    (
        this->lookup("twoFluids")
    ),
    normalMotionDir_
    (
        this->lookup("normalMotionDir")
    ),
    motionDir_(0, 0, 0),
    cleanInterface_
    (
        this->lookup("cleanInterface")
    ),
    aPatchID_(-1),
    bPatchID_(-1),
    muFluidA_
    (
        transportModel.rho1()
      * dimensionedScalar(transportModel.nuModel1().viscosityProperties().lookup("nu"))
    ),
    muFluidB_
    (
        transportModel.rho2()
      * dimensionedScalar(transportModel.nuModel2().viscosityProperties().lookup("nu"))
    ),
    rhoFluidA_
    (
        transportModel.rho1()
    ),
    rhoFluidB_
    (
        transportModel.rho2()
    ),
    cleanInterfaceSurfTension_
    (
        dimensionedScalar(transportModel.lookup("sigma"))
    ),
    fixedTrackedSurfacePatches_
    (
        this->lookup("fixedTrackedSurfacePatches")
    ),
    pointNormalsCorrectionPatches_
    (
        this->lookup("pointNormalsCorrectionPatches")
    ),
    nTrackedSurfCorr_
    (
        readInt(this->lookup("nTrackedSurfaceCorrectors"))
    ),
    smoothing_(false),
    noCleanTangentialSurfaceTensionCorrection_
    (
        this->lookupOrDefault<Switch>("noCleanTangentialSurfaceTensionCorrection", false)
    ),
    interpolatorABPtr_(NULL),
    interpolatorBAPtr_(NULL),
    controlPointsPtr_(NULL),
    motionPointsMaskPtr_(NULL),
    pointsDisplacementDirPtr_(NULL),
    facesDisplacementDirPtr_(NULL),
    totalDisplacementPtr_(NULL),
    aMeshPtr_(NULL),
    UsPtr_(NULL),
    phisPtr_(NULL),
    surfactConcPtr_(NULL),
    surfaceTensionPtr_(NULL),
    surfactantPtr_(NULL),
    fluidIndicatorPtr_(NULL),
    muEffFluidAvalPtr_(NULL),
    muEffFluidBvalPtr_(NULL)
{
    //Read motion direction
    if (!normalMotionDir_)
    {
        motionDir_ = vector(this->lookup("motionDir"));
        motionDir_ /= mag(motionDir_) + SMALL;
    }

    // Set point normal correction patches
    boolList& correction = aMesh().correctPatchPointNormals();

    forAll(pointNormalsCorrectionPatches_, patchI)
    {
        word patchName = pointNormalsCorrectionPatches_[patchI];

        label patchID = aMesh().boundary().findPatchID(patchName);

        if(patchID == -1)
        {
            FatalErrorIn
            (
                "trackedSurface::trackedSurface(...)"
            )   << "Patch name for point normals correction does not exist"
                << abort(FatalError);
        }

        correction[patchID] = true;
    }

    // Clear geometry
    aMesh().movePoints();


    // Detect the surface patch
    forAll (mesh().boundary(), patchI)
    {
        if(mesh().boundary()[patchI].name() == "trackedSurface")
        {
            aPatchID_ = patchI;

            Info<< "Found surface patch. ID: " << aPatchID_
                << endl;
        }
    }

    if(aPatchID() == -1)
    {
        FatalErrorIn("trackedSurface::trackedSurface(...)")
            << "Surface patch not defined.  Please make sure that "
                << " the surface patches is named as trackedSurface"
                << abort(FatalError);
    }


    // Detect the surface shadow patch
    if (twoFluids())
    {
        forAll (mesh().boundary(), patchI)
        {
            if(mesh().boundary()[patchI].name() == "trackedSurfaceShadow")
            {
                bPatchID_ = patchI;

                Info<< "Found surface shadow patch. ID: "
                    << bPatchID_ << endl;
            }
        }

        if(bPatchID() == -1)
        {
            FatalErrorIn("trackedSurface::trackedSurface(...)")
                << "Surface shadow patch not defined. "
                    << "Please make sure that the surface shadow patch "
                    << "is named as trackedSurfaceShadow."
                    << abort(FatalError);
        }
    }


    // Mark surface boundary points
    // which belonge to processor patches
    forAll(aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type()
         == processorFaPatch::typeName
        )
        {
            const labelList& patchPoints =
                aMesh().boundary()[patchI].pointLabels();

            forAll(patchPoints, pointI)
            {
                motionPointsMask()[patchPoints[pointI]] = -1;
            }
        }
    }


    // Mark fixed surface boundary points
    forAll(fixedTrackedSurfacePatches_, patchI)
    {
        label fixedPatchID =
            aMesh().boundary().findPatchID
            (
                fixedTrackedSurfacePatches_[patchI]
            );

        if(fixedPatchID == -1)
        {
            FatalErrorIn("trackedSurface::trackedSurface(...)")
                << "Wrong faPatch name in the fixedTrackedSurfacePatches list"
                    << " defined in the trackedSurfaceProperties dictionary"
                    << abort(FatalError);
        }

        const labelList& patchPoints =
            aMesh().boundary()[fixedPatchID].pointLabels();

        forAll(patchPoints, pointI)
        {
            motionPointsMask()[patchPoints[pointI]] = 0;
        }
    }


    // Mark surface boundary point
    // at the axis of 2-D axisymmetic cases
    forAll(aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type()
         == wedgeFaPatch::typeName
        )
        {
            const wedgeFaPatch& wedgePatch =
                refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

            if(wedgePatch.axisPoint() > -1)
            {
                motionPointsMask()[wedgePatch.axisPoint()] = 0;

                Info << "Axis point: "
                    << wedgePatch.axisPoint()
                    << "vector: "
                    << aMesh().points()[wedgePatch.axisPoint()] << endl;
            }
        }
    }


    // Read surface points total displacement if present
    readTotalDisplacement();


    // Read control points positions if present
    controlPoints();


    // Check if smoothing switch is set
    if (this->found("smoothing"))
    {
        smoothing_ = Switch(this->lookup("smoothing"));
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

trackedSurface::~trackedSurface()
{
    clearOut();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void trackedSurface::updateDisplacementDirections()
{
    if(normalMotionDir())
    {
        // Update point displacement correction
        pointsDisplacementDir() = aMesh().pointAreaNormals();

        // Correcte point displacement direction
        // at the "centerline" symmetryPlane which represents the axis
        // of an axisymmetric case
        forAll(aMesh().boundary(), patchI)
        {
            if(aMesh().boundary()[patchI].type() == wedgeFaPatch::typeName)
            {
                const wedgeFaPatch& wedgePatch =
                    refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

                vector axis = wedgePatch.axis();

                label centerLinePatchID =
                    aMesh().boundary().findPatchID("centerline");

                if(centerLinePatchID != -1)
                {
                    const labelList& pointLabels =
                        aMesh().boundary()[centerLinePatchID].pointLabels();

                    forAll(pointLabels, pointI)
                    {
                        vector dir =
                            pointsDisplacementDir()[pointLabels[pointI]];

                        dir = (dir&axis)*axis;
                        dir /= mag(dir);

                        pointsDisplacementDir()[pointLabels[pointI]] = dir;
                    }
                }
                else
                {
                    Info << "Warning: centerline polyPatch does not exist. "
                        << "Surface points displacement directions "
                        << "will not be corrected at the axis (centerline)"
                        << endl;
                }

                break;
            }
        }

        // Update face displacement direction
        facesDisplacementDir() =
            aMesh().faceAreaNormals().internalField();

        // Correction of control points postion
        const vectorField& Cf = aMesh().areaCentres().internalField();

        controlPoints() =
            facesDisplacementDir()
           *(facesDisplacementDir()&(controlPoints() - Cf))
          + Cf;
    }
}


bool trackedSurface::predictPoints()
{
    // Smooth interface

    if (smoothing_)
    {
        controlPoints() = aMesh().areaCentres().internalField();
        movePoints(scalarField(controlPoints().size(), 0));
        movePoints(-fvc::meshPhi(U())().boundaryField()[aPatchID()]);
    }

    for
    (
        int trackedSurfCorr=0;
        trackedSurfCorr<nTrackedSurfCorr_;
        trackedSurfCorr++
    )
    {
        movePoints(phi_.boundaryField()[aPatchID()]);
    }

    return true;
}


bool trackedSurface::correctPoints()
{
    for
    (
        int trackedSurfCorr=0;
        trackedSurfCorr<nTrackedSurfCorr_;
        trackedSurfCorr++
    )
    {
        movePoints(phi_.boundaryField()[aPatchID()]);
    }

    return true;
}


bool trackedSurface::movePoints(const scalarField& interfacePhi)
{
    pointField newMeshPoints = mesh().points();

    scalarField sweptVolCorr =
        interfacePhi
      - fvc::meshPhi(U())().boundaryField()[aPatchID()];

    word ddtScheme
    (
        mesh().schemesDict().ddtScheme
        (
            "ddt(" + U().name() + ")"
        )
    );

    if
    (
        ddtScheme
     == fv::CrankNicolsonDdtScheme<vector>::typeName
    )
    {
        sweptVolCorr *= (1.0/2.0)*DB().deltaT().value();
    }
    else if
    (
        ddtScheme
     == fv::EulerDdtScheme<vector>::typeName
    )
    {
        sweptVolCorr *= DB().deltaT().value();
    }
    else if
    (
        ddtScheme
     == fv::backwardDdtScheme<vector>::typeName
    )
    {
        if (DB().timeIndex() == 1)
        {
            sweptVolCorr *= DB().deltaT().value();
        }
        else
        {
            sweptVolCorr *= (2.0/3.0)*DB().deltaT().value();
        }
    }
    else
    {
        FatalErrorIn("trackedSurface::movePoints()")
            << "Unsupported temporal differencing scheme : "
                << ddtScheme
                << abort(FatalError);
    }

    const scalarField& Sf = aMesh().S();
    const vectorField& Nf = aMesh().faceAreaNormals().internalField();

    scalarField deltaH =
        sweptVolCorr/(Sf*(Nf & facesDisplacementDir()));

    pointField displacement = pointDisplacement(deltaH);


    // Move only surface points

    const labelList& meshPointsA =
        mesh().boundaryMesh()[aPatchID()].meshPoints();

    forAll (displacement, pointI)
    {
        newMeshPoints[meshPointsA[pointI]] += displacement[pointI];
    }

    if(twoFluids_)
    {
        const labelList& meshPointsB =
            mesh().boundaryMesh()[bPatchID_].meshPoints();

        pointField displacementB =
            interpolatorAB().pointInterpolate
            (
                displacement
            );

        forAll (displacementB, pointI)
        {
            newMeshPoints[meshPointsB[pointI]] += displacementB[pointI];
        }
    }

    // Update total displacement field

    if(totalDisplacementPtr_ && (curTimeIndex_ < DB().timeIndex()))
    {
        FatalErrorIn("trackedSurface::movePoints()")
            << "Total displacement of surface points "
                << "from previous time step is not absorbed by the mesh."
                << abort(FatalError);
    }
    else if (curTimeIndex_ < DB().timeIndex())
    {
        totalDisplacement() = displacement;

        curTimeIndex_ = DB().timeIndex();
    }
    else
    {
        totalDisplacement() += displacement;
    }

    twoDPointCorrector twoDPointCorr(mesh());

    twoDPointCorr.correctPoints(newMeshPoints);

    mesh().movePoints(newMeshPoints);

    // faMesh motion is done automatically, using meshObject
    // HJ, 8/Aug/2011
//     aMesh().movePoints(mesh().points());


    // Move correctedFvPatchField fvSubMeshes

    forAll(U().boundaryField(), patchI)
    {
        if
        (
            (
                U().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<vector>::typeName
            )
            ||
            (
                U().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<vector>::typeName
            )
            ||
            (
                U().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<vector>::typeName
            )
        )
        {
            correctedFvPatchField<vector>& pU =
                refCast<correctedFvPatchField<vector> >
                (
                    U().boundaryField()[patchI]
                );

            pU.movePatchSubMesh();
        }
    }

    forAll(p().boundaryField(), patchI)
    {
        if
        (
            (
                p().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<scalar>::typeName
            )
            ||
            (
                p().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<scalar>::typeName
            )
            ||
            (
                p().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<scalar>::typeName
            )
        )
        {
            correctedFvPatchField<scalar>& pP =
                refCast<correctedFvPatchField<scalar> >
                (
                    p().boundaryField()[patchI]
                );

            pP.movePatchSubMesh();
        }
    }

    return true;
}


bool trackedSurface::moveMeshPointsForOldTrackedSurfDisplacement()
{
    if(totalDisplacementPtr_)
    {
        pointField newPoints = mesh().points();

        const labelList& meshPointsA =
            mesh().boundaryMesh()[aPatchID()].meshPoints();

        forAll (totalDisplacement(), pointI)
        {
            newPoints[meshPointsA[pointI]] -= totalDisplacement()[pointI];
        }


        // Check mesh motion solver type
        bool feMotionSolver =
            mesh().objectRegistry::foundObject<tetPointVectorField>
            (
                "motionU"
            );
        bool fvMotionSolver =
            mesh().objectRegistry::foundObject<pointVectorField>
            (
                "pointMotionU"
            );

        if (feMotionSolver)
        {
            tetPointVectorField& motionU =
                const_cast<tetPointVectorField&>
                (
                    mesh().objectRegistry::
                    lookupObject<tetPointVectorField>
                    (
                        "motionU"
                    )
                );

            fixedValueTetPolyPatchVectorField& motionUaPatch =
                refCast<fixedValueTetPolyPatchVectorField>
                (
                    motionU.boundaryField()[aPatchID()]
                );

            tetPolyPatchInterpolation tppiAPatch
            (
                refCast<const faceTetPolyPatch>
                (
                    motionUaPatch.patch()
                )
            );

            motionUaPatch ==
                tppiAPatch.pointToPointInterpolate
                (
                    totalDisplacement()/DB().deltaT().value()
                );

            if(twoFluids_)
            {
                const labelList& meshPointsB =
                    mesh().boundaryMesh()[bPatchID()].meshPoints();

                pointField totDisplacementB =
                    interpolatorAB().pointInterpolate
                    (
                        totalDisplacement()
                    );

                forAll (totDisplacementB, pointI)
                {
                    newPoints[meshPointsB[pointI]] -=
                        totDisplacementB[pointI];
                }

                fixedValueTetPolyPatchVectorField& motionUbPatch =
                    refCast<fixedValueTetPolyPatchVectorField>
                    (
                        motionU.boundaryField()[bPatchID()]
                    );

                tetPolyPatchInterpolation tppiBPatch
                (
                    refCast<const faceTetPolyPatch>(motionUbPatch.patch())
                );

                motionUbPatch ==
                    tppiBPatch.pointToPointInterpolate
                    (
                        totDisplacementB/DB().deltaT().value()
                    );
            }
        }
        else if (fvMotionSolver)
        {
            pointVectorField& motionU =
                const_cast<pointVectorField&>
                (
                    mesh().objectRegistry::
                    lookupObject<pointVectorField>
                    (
                        "pointMotionU"
                    )
                );

            fixedValuePointPatchVectorField& motionUaPatch =
                refCast<fixedValuePointPatchVectorField>
                (
                    motionU.boundaryField()[aPatchID()]
                );

            motionUaPatch ==
                totalDisplacement()/DB().deltaT().value();

            if(twoFluids_)
            {
                const labelList& meshPointsB =
                    mesh().boundaryMesh()[bPatchID()].meshPoints();

                pointField totDisplacementB =
                    interpolatorAB().pointInterpolate
                    (
                        totalDisplacement()
                    );

                forAll (totDisplacementB, pointI)
                {
                    newPoints[meshPointsB[pointI]] -=
                        totDisplacementB[pointI];
                }

                fixedValuePointPatchVectorField& motionUbPatch =
                    refCast<fixedValuePointPatchVectorField>
                    (
                        motionU.boundaryField()[bPatchID()]
                    );

                motionUbPatch ==
                    totDisplacementB/DB().deltaT().value();
            }
        }

        twoDPointCorrector twoDPointCorr(mesh());

        twoDPointCorr.correctPoints(newPoints);

        mesh().movePoints(newPoints);

        deleteDemandDrivenData(totalDisplacementPtr_);

        mesh().update();

        // faMesh motion is done automatically, using meshObject
        // HJ, 8/Aug/2011
//         aMesh().movePoints(mesh().points());

        // Move correctedFvPatchField fvSubMeshes

        forAll(U().boundaryField(), patchI)
        {
            if
            (
                (
                    U().boundaryField()[patchI].type()
                 == fixedGradientCorrectedFvPatchField<vector>::typeName
                )
                ||
                (
                    U().boundaryField()[patchI].type()
                 == fixedValueCorrectedFvPatchField<vector>::typeName
                )
                ||
                (
                    U().boundaryField()[patchI].type()
                 == zeroGradientCorrectedFvPatchField<vector>::typeName
                )
            )
            {
                correctedFvPatchField<vector>& aU =
                    refCast<correctedFvPatchField<vector> >
                    (
                        U().boundaryField()[patchI]
                    );

                aU.movePatchSubMesh();
            }
        }

        forAll(p().boundaryField(), patchI)
        {
            if
            (
                (
                    p().boundaryField()[patchI].type()
                 == fixedGradientCorrectedFvPatchField<scalar>::typeName
                )
                ||
                (
                    p().boundaryField()[patchI].type()
                 == fixedValueCorrectedFvPatchField<scalar>::typeName
                )
                ||
                (
                    p().boundaryField()[patchI].type()
                 == zeroGradientCorrectedFvPatchField<scalar>::typeName
                )
            )
            {
                correctedFvPatchField<scalar>& aP =
                    refCast<correctedFvPatchField<scalar> >
                    (
                        p().boundaryField()[patchI]
                    );

                aP.movePatchSubMesh();
            }
        }
    }

    return true;
}


bool trackedSurface::moveMeshPoints()
{
        scalarField sweptVolCorr =
            phi_.boundaryField()[aPatchID()]
          - fvc::meshPhi(U())().boundaryField()[aPatchID()];

        word ddtScheme
        (
            mesh().schemesDict().ddtScheme
            (
                "ddt(" + U().name() + ")"
            )
        );

        if
        (
            ddtScheme
         == fv::CrankNicolsonDdtScheme<vector>::typeName
        )
        {
            sweptVolCorr *= (1.0/2.0)*DB().deltaT().value();
        }
        else if
        (
            ddtScheme
         == fv::EulerDdtScheme<vector>::typeName
        )
        {
            sweptVolCorr *= DB().deltaT().value();
        }
        else if
        (
            ddtScheme
         == fv::backwardDdtScheme<vector>::typeName
        )
        {
            sweptVolCorr *= (2.0/3.0)*DB().deltaT().value();
        }
        else
        {
            FatalErrorIn("trackedSurface::movePoints()")
                << "Unsupported temporal differencing scheme : "
                << ddtScheme
                << abort(FatalError);
        }


        const scalarField& Sf = aMesh().S();
        const vectorField& Nf = aMesh().faceAreaNormals().internalField();

        scalarField deltaH =
            sweptVolCorr/(Sf*(Nf & facesDisplacementDir()));


        pointField displacement = pointDisplacement(deltaH);


        //-- Set mesh motion boundary conditions

        tetPointVectorField& motionU =
            const_cast<tetPointVectorField&>
            (
                mesh().objectRegistry::
                lookupObject<tetPointVectorField>
                (
                    "motionU"
                )
            );

        fixedValueTetPolyPatchVectorField& motionUaPatch =
            refCast<fixedValueTetPolyPatchVectorField>
            (
                motionU.boundaryField()[aPatchID()]
            );

        tetPolyPatchInterpolation tppiAPatch
        (
            refCast<const faceTetPolyPatch>
            (
                motionUaPatch.patch()
            )
        );

        motionUaPatch ==
            tppiAPatch.pointToPointInterpolate
            (
                displacement/DB().deltaT().value()
            );

        if (twoFluids())
        {
            fixedValueTetPolyPatchVectorField& motionUbPatch =
                refCast<fixedValueTetPolyPatchVectorField>
                (
                    motionU.boundaryField()[bPatchID()]
                );

            tetPolyPatchInterpolation tppiBPatch
            (
                refCast<const faceTetPolyPatch>(motionUbPatch.patch())
            );

            motionUbPatch ==
                tppiBPatch.pointToPointInterpolate
                (
                    interpolatorAB().pointInterpolate
                    (
                        displacement/DB().deltaT().value()
                    )
                );
        }

        mesh().update();

        // faMesh motion is done automatically, using meshObject
        // HJ, 8/Aug/2011
//         aMesh().movePoints(mesh().points());


        // Move correctedFvPatchField fvSubMeshes

        forAll(U().boundaryField(), patchI)
        {
            if
            (
                (
                    U().boundaryField()[patchI].type()
                 == fixedGradientCorrectedFvPatchField<vector>::typeName
                )
                ||
                (
                    U().boundaryField()[patchI].type()
                 == fixedValueCorrectedFvPatchField<vector>::typeName
                )
                ||
                (
                    U().boundaryField()[patchI].type()
                 == zeroGradientCorrectedFvPatchField<vector>::typeName
                )
            )
            {
                correctedFvPatchField<vector>& aU =
                    refCast<correctedFvPatchField<vector> >
                    (
                        U().boundaryField()[patchI]
                    );

                aU.movePatchSubMesh();
            }
        }

        forAll(p().boundaryField(), patchI)
        {
            if
            (
                (
                    p().boundaryField()[patchI].type()
                 == fixedGradientCorrectedFvPatchField<scalar>::typeName
                )
                ||
                (
                    p().boundaryField()[patchI].type()
                 == fixedValueCorrectedFvPatchField<scalar>::typeName
                )
                ||
                (
                    p().boundaryField()[patchI].type()
                 == zeroGradientCorrectedFvPatchField<scalar>::typeName
                )
            )
            {
                correctedFvPatchField<scalar>& aP =
                    refCast<correctedFvPatchField<scalar> >
                    (
                        p().boundaryField()[patchI]
                    );

                aP.movePatchSubMesh();
            }
        }

    return true;
}


void trackedSurface::updateBoundaryConditions()
{
    updateMuEff();
    updateVelocity();
    updateSurfactantConcentration();
    updatePressure();
}


void trackedSurface::updateMuEff()
{
    const volScalarField& nuEff = turbulence()->nuEff();

    muEffFluidAval() =
        nuEff.boundaryField()[aPatchID()] * rhoFluidA().value();

    if(twoFluids())
    {
        muEffFluidBval() =
            nuEff.boundaryField()[bPatchID()] * rhoFluidB().value();
    }
}


void trackedSurface::updateVelocity()
{
    if(twoFluids())
    {
        vectorField nA = mesh().boundary()[aPatchID()].nf();

        vectorField nB = mesh().boundary()[bPatchID()].nf();

        scalarField DnB = interpolatorBA().faceInterpolate
        (
            mesh().boundary()[bPatchID()].deltaCoeffs()
        );

        scalarField DnA = mesh().boundary()[aPatchID()].deltaCoeffs();


        vectorField UtPA =
            U().boundaryField()[aPatchID()].patchInternalField();

        if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientCorrectedFvPatchField<vector>::typeName
        )
        {
            fixedGradientCorrectedFvPatchField<vector>& aU =
                refCast<fixedGradientCorrectedFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );

            UtPA += aU.corrVecGrad();
        }

        UtPA -= nA*(nA & UtPA);


        vectorField UtPB = interpolatorBA().faceInterpolate
        (
            U().boundaryField()[bPatchID()].patchInternalField()
        );

        if
        (
            U().boundaryField()[bPatchID()].type()
         == fixedValueCorrectedFvPatchField<vector>::typeName
        )
        {
            fixedValueCorrectedFvPatchField<vector>& bU =
                refCast<fixedValueCorrectedFvPatchField<vector> >
                (
                    U().boundaryField()[bPatchID()]
                );

            UtPB += interpolatorBA().faceInterpolate(bU.corrVecGrad());
        }

        UtPB -= nA*(nA & UtPB);

        vectorField UtFs = muEffFluidAval()*DnA*UtPA
          + muEffFluidBval()*DnB*UtPB;

        vectorField UnFs =
            nA*phi_.boundaryField()[aPatchID()]
           /mesh().boundary()[aPatchID()].magSf();

        Us().internalField() += UnFs - nA*(nA&Us().internalField());
        correctUsBoundaryConditions();

        UtFs -= (muEffFluidAval() - muEffFluidBval())*
            (fac::grad(Us())&aMesh().faceAreaNormals())().internalField();


        vectorField tangentialSurfaceTensionForce(nA.size(), vector::zero);

        tangentialSurfaceTensionForce =
            surfaceTensionGrad()().internalField();

        UtFs += tangentialSurfaceTensionForce;

        UtFs /= muEffFluidAval()*DnA + muEffFluidBval()*DnB + VSMALL;

        Us().internalField() = UnFs + UtFs;
        correctUsBoundaryConditions();

        // Store old-time velocity field U()
        U().oldTime();

        U().boundaryField()[bPatchID()] ==
            interpolatorAB().faceInterpolate(UtFs)
          + nB*fvc::meshPhi(U())().boundaryField()[bPatchID()]/
            mesh().boundary()[bPatchID()].magSf();

        if
        (
            p().boundaryField()[bPatchID()].type()
         == fixedGradientFvPatchField<scalar>::typeName
        )
        {
            fixedGradientFvPatchField<scalar>& pB =
                refCast<fixedGradientFvPatchField<scalar> >
                (
                    p().boundaryField()[bPatchID()]
                );

            pB.gradient() =
               - rhoFluidB().value()
                *(
                     nB&fvc::ddt(U())().boundaryField()[bPatchID()]
                 );
        }


        // Update fixedGradient boundary condition on patch A

        vectorField nGradU =
            muEffFluidBval()*(UtPB - UtFs)*DnA
          + tangentialSurfaceTensionForce
          - muEffFluidAval()*nA*fac::div(Us())().internalField()
          + (muEffFluidBval() - muEffFluidAval())
           *(fac::grad(Us())().internalField()&nA);

        nGradU /= muEffFluidAval() + VSMALL;

        if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientCorrectedFvPatchField<vector>::typeName
        )
        {
            fixedGradientCorrectedFvPatchField<vector>& aU =
                refCast<fixedGradientCorrectedFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );

            aU.gradient() = nGradU;
        }
        else if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientFvPatchField<vector>::typeName
        )
        {
            fixedGradientFvPatchField<vector>& aU =
                refCast<fixedGradientFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );

            aU.gradient() = nGradU;
        }
        else
        {
            FatalErrorIn("trackedSurface::updateVelocity()")
                << "Bounary condition on " << U().name()
                    <<  " for trackedSurface patch is "
                    << U().boundaryField()[aPatchID()].type()
                    << ", instead "
                    << fixedGradientCorrectedFvPatchField<vector>::typeName
                    << " or "
                    << fixedGradientFvPatchField<vector>::typeName
                    << abort(FatalError);
        }
    }
    else
    {
        vectorField nA = aMesh().faceAreaNormals().internalField();

        vectorField UnFs =
            nA*phi_.boundaryField()[aPatchID()]
           /mesh().boundary()[aPatchID()].magSf();

        // Correct normal component of surface velocity
        Us().internalField() += UnFs - nA*(nA&Us().internalField());
        correctUsBoundaryConditions();

        vectorField tangentialSurfaceTensionForce(nA.size(), vector::zero);

        tangentialSurfaceTensionForce =
            surfaceTensionGrad()().internalField();

        vectorField tnGradU =
            tangentialSurfaceTensionForce/(muEffFluidAval() + VSMALL)
          - (fac::grad(Us())&aMesh().faceAreaNormals())().internalField();

        vectorField UtPA =
            U().boundaryField()[aPatchID()].patchInternalField();
        UtPA -= nA*(nA & UtPA);

        scalarField DnA = mesh().boundary()[aPatchID()].deltaCoeffs();

        vectorField UtFs = UtPA + tnGradU/DnA;

        Us().internalField() = UtFs + UnFs;
        correctUsBoundaryConditions();

        vectorField nGradU =
            tangentialSurfaceTensionForce/(muEffFluidAval() + VSMALL)
          - nA*fac::div(Us())().internalField()
          - (fac::grad(Us())().internalField()&nA);

        if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientCorrectedFvPatchField<vector>::typeName
        )
        {
            fixedGradientCorrectedFvPatchField<vector>& aU =
                refCast<fixedGradientCorrectedFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );

            aU.gradient() = nGradU;
        }
        else if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientFvPatchField<vector>::typeName
        )
        {
            fixedGradientFvPatchField<vector>& aU =
                refCast<fixedGradientFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );

            aU.gradient() = nGradU;
        }
        else
        {
            FatalErrorIn("trackedSurface::updateVelocity()")
                << "Bounary condition on " << U().name()
                    <<  " for trackedSurface patch is "
                    << U().boundaryField()[aPatchID()].type()
                    << ", instead "
                    << fixedGradientCorrectedFvPatchField<vector>::typeName
                    << " or "
                    << fixedGradientFvPatchField<vector>::typeName
                    << abort(FatalError);
        }
    }
}


void trackedSurface::updatePressure()
{
    // Correct pressure boundary condition at the surface

    vectorField nA = mesh().boundary()[aPatchID()].nf();

    if(twoFluids())
    {
        scalarField pA =
            interpolatorBA().faceInterpolate
            (
                p().boundaryField()[bPatchID()]
            );

        const scalarField& K = aMesh().faceCurvatures().internalField();

        Info << "Surface curvature: min = " << gMin(K)
            << ", max = " << gMax(K)
            << ", average = " << gAverage(K) << endl << flush;

        if(cleanInterface())
        {
//             pA -= cleanInterfaceSurfTension().value()*(K - gAverage(K));
            pA -= cleanInterfaceSurfTension().value()*K;
        }
        else
        {
            scalarField surfTensionK =
                surfaceTension().internalField()*K;

            pA -= surfTensionK - gAverage(surfTensionK);
        }

        pA -= 2.0*(muEffFluidAval() - muEffFluidBval())
            *fac::div(Us())().internalField();

        pA += pExt().boundaryField()[aPatchID()];

//         vector R0 = gAverage(mesh().C().boundaryField()[aPatchID()]);
        vector R0 = vector::zero;

        pA -= (rhoFluidA().value() - rhoFluidB().value())*
            (
                g_.value()
              & (
                    mesh().C().boundaryField()[aPatchID()]
                  - R0
                )
            );

        p().boundaryField()[aPatchID()] == pA;
    }
    else
    {
//         vector R0 = gAverage(mesh().C().boundaryField()[aPatchID()]);
        vector R0 = vector::zero;

        scalarField pA =
          - rhoFluidA().value()*
            (
                g_.value()
              & (
                    mesh().C().boundaryField()[aPatchID()]
                  - R0
                )
            );

        const scalarField& K = aMesh().faceCurvatures().internalField();

        Info << "Surface curvature: min = " << gMin(K)
            << ", max = " << gMax(K) << ", average = " << gAverage(K)
            << endl;

        if(cleanInterface())
        {
//             pA -= cleanInterfaceSurfTension().value()*(K - gAverage(K));
            pA -= cleanInterfaceSurfTension().value()*K;
        }
        else
        {
            scalarField surfTensionK =
                surfaceTension().internalField()*K;

            pA -= surfTensionK - gAverage(surfTensionK);
        }

        pA -= 2.0*muEffFluidAval()*fac::div(Us())().internalField();

        pA += pExt().boundaryField()[aPatchID()];

        p().boundaryField()[aPatchID()] == pA;
    }


    // Set modified pressure at patches with fixed apsolute
    // pressure

//     vector R0 = gAverage(mesh().C().boundaryField()[aPatchID()]);
    vector R0 = vector::zero;

    for (int patchI=0; patchI < p().boundaryField().size(); patchI++)
    {
        if
        (
            p().boundaryField()[patchI].type()
         == fixedValueFvPatchScalarField::typeName
        )
        {
            if (patchI != aPatchID())
            {
                p().boundaryField()[patchI] ==
                    pExt().boundaryField()[patchI];
                  - rho().boundaryField()[patchI]
                   *(g_.value()&(mesh().C().boundaryField()[patchI] - R0));
            }
        }
    }
}


void trackedSurface::updateSurfaceFlux()
{
    Phis() = fac::interpolate(Us()) & aMesh().Le();
}


void trackedSurface::updateSurfactantConcentration()
{
    if(!cleanInterface())
    {
        Info << "Correct surfactant concentration" << endl << flush;

        updateSurfaceFlux();

        // Crate and solve the surfactanta transport equation
        faScalarMatrix CsEqn
        (
            fam::ddt(surfactantConcentration())
          + fam::div(Phis(), surfactantConcentration())
          - fam::laplacian
            (
                surfactant().surfactDiffusion(),
                surfactantConcentration()
            )
        );


        if(surfactant().soluble())
        {
            const scalarField& C =
                mesh().boundary()[aPatchID()]
               .lookupPatchField<volScalarField, scalar>("C");

            areaScalarField Cb
            (
                IOobject
                (
                    "Cb",
                    DB().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                aMesh(),
                dimensioned<scalar>("Cb", dimMoles/dimVolume, 0),
                zeroGradientFaPatchScalarField::typeName
            );

            Cb.internalField() = C;
            Cb.correctBoundaryConditions();

            CsEqn +=
                fam::Sp
                (
                    surfactant().surfactAdsorptionCoeff()*Cb
                  + surfactant().surfactAdsorptionCoeff()
                   *surfactant().surfactDesorptionCoeff(),
                    surfactantConcentration()
                )
              - surfactant().surfactAdsorptionCoeff()
               *Cb*surfactant().surfactSaturatedConc();
        }

        CsEqn.solve();

        Info << "Correct surface tension" << endl;

        surfaceTension() =
            cleanInterfaceSurfTension()
          + surfactant().surfactR()
           *surfactant().surfactT()
           *surfactant().surfactSaturatedConc()
           *log(1.0 - surfactantConcentration()
           /surfactant().surfactSaturatedConc());

        if(neg(min(surfaceTension().internalField())))
        {
            FatalErrorIn
            (
                "void trackedSurface::correctSurfactantConcentration()"
            )   << "Surface tension is negative"
                    << abort(FatalError);
        }
    }
}


void trackedSurface::correctUsBoundaryConditions()
{
    forAll(Us().boundaryField(), patchI)
    {
        if
        (
            Us().boundaryField()[patchI].type()
         == calculatedFaPatchVectorField::typeName
        )
        {
            vectorField& pUs = Us().boundaryField()[patchI];

            pUs = Us().boundaryField()[patchI].patchInternalField();

            label ngbPolyPatchID =
                aMesh().boundary()[patchI].ngbPolyPatchIndex();

            if(ngbPolyPatchID != -1)
            {
                if
                (
                    (
                        isA<slipFvPatchVectorField>
                        (
                            U().boundaryField()[ngbPolyPatchID]
                        )
                    )
                 ||
                    (
                        isA<symmetryFvPatchVectorField>
                        (
                            U().boundaryField()[ngbPolyPatchID]
                        )
                    )
                )
                {
                    vectorField N =
                        aMesh().boundary()[patchI].ngbPolyPatchFaceNormals();

                    pUs -= N*(N&pUs);
                }
            }
        }
    }

    Us().correctBoundaryConditions();
}


vector trackedSurface::totalPressureForce() const
{
    const scalarField& S = aMesh().S();

    const vectorField& n = aMesh().faceAreaNormals().internalField();

    const scalarField& P = p().boundaryField()[aPatchID()];

    vectorField pressureForces = S*P*n;

    return gSum(pressureForces);
}


vector trackedSurface::totalViscousForce() const
{
    const scalarField& S = aMesh().S();
    const vectorField& n = aMesh().faceAreaNormals().internalField();

    vectorField nGradU =
        U().boundaryField()[aPatchID()].snGrad();

    vectorField viscousForces =
      - muEffFluidAval()*S
       *(
            nGradU
          + (fac::grad(Us())().internalField()&n)
          - (n*fac::div(Us())().internalField())
        );

    return gSum(viscousForces);
}


vector trackedSurface::totalSurfaceTensionForce() const
{
    const scalarField& S = aMesh().S();

    const vectorField& n = aMesh().faceAreaNormals().internalField();

    const scalarField& K = aMesh().faceCurvatures().internalField();

    vectorField surfTensionForces(n.size(), vector::zero);

    if(cleanInterface())
    {
        surfTensionForces =
            S*cleanInterfaceSurfTension().value()
           *fac::edgeIntegrate
            (
                aMesh().Le()*aMesh().edgeLengthCorrection()
            )().internalField();
    }
    else
    {
        surfTensionForces *= surfaceTension().internalField()*K;
    }

    return gSum(surfTensionForces);
}


void trackedSurface::initializeControlPointsPosition()
{
    scalarField deltaH = scalarField(aMesh().nFaces(), 0.0);

    pointField displacement = pointDisplacement(deltaH);

    const faceList& faces = aMesh().faces();
    const pointField& points = aMesh().points();


    pointField newPoints = points + displacement;

    scalarField sweptVol(faces.size(), 0.0);

    forAll(faces, faceI)
    {
        sweptVol[faceI] = -faces[faceI].sweptVol(points, newPoints);
    }

    vectorField faceArea(faces.size(), vector::zero);

    forAll (faceArea, faceI)
    {
        faceArea[faceI] = faces[faceI].normal(newPoints);
    }

    forAll(deltaH, faceI)
    {
        deltaH[faceI] = sweptVol[faceI]/
            (faceArea[faceI] & facesDisplacementDir()[faceI]);
    }

    displacement = pointDisplacement(deltaH);
}


scalar trackedSurface::maxCourantNumber()
{
    scalar CoNum = 0;

    if(cleanInterface())
    {
        const scalarField& dE =aMesh().lPN();

        CoNum = gMax
        (
            DB().deltaT().value()/
            sqrt
            (
                rhoFluidA().value()*dE*dE*dE/
                2.0/M_PI/(cleanInterfaceSurfTension().value() + SMALL)
            )
        );
    }
    else
    {
        scalarField sigmaE =
            linearEdgeInterpolate(surfaceTension())().internalField()
          + SMALL;

        const scalarField& dE =aMesh().lPN();

        CoNum = gMax
        (
            DB().deltaT().value()/
            sqrt
            (
                rhoFluidA().value()*dE*dE*dE/
                2.0/M_PI/sigmaE
            )
        );
    }

    return CoNum;
}


void trackedSurface::updateProperties()
{
    muFluidA_ = transport().rho1()
      * dimensionedScalar
        (
            transport().nuModel1().viscosityProperties().lookup("nu")
        );

    muFluidB_ = transport().rho2()
      * dimensionedScalar
        (
            transport().nuModel2().viscosityProperties().lookup("nu")
        );

    rhoFluidA_ = transport().rho1();

    rhoFluidB_ = transport().rho2();

    cleanInterfaceSurfTension_ =
        dimensionedScalar(transport().lookup("sigma"));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
