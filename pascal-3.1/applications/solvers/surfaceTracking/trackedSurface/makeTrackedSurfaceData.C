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

\*---------------------------------------------------------------------------*/

#include "trackedSurface.H"
#include "primitivePatchInterpolation.H"
#include "wedgeFaPatch.H"
#include "wallFvPatch.H"
#include "wedgeFaPatchFields.H"
#include "slipFaPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void trackedSurface::makeInterpolators()
{
    if (debug)
    {
        Info<< "trackedSurface::makeInterpolators() : "
            << "Making patch to patch interpolator"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if
    (
        interpolatorBAPtr_ ||
        interpolatorABPtr_
    )
    {
        FatalErrorIn("trackedSurface::makeInterpolators()")
            << "patch to patch interpolators already exists"
                << abort(FatalError);
    }


    if(aPatchID() == -1)
    {
        FatalErrorIn("trackedSurface::makeInterpolators()")
            << "Surface patch A not defined."
            << abort(FatalError);
    }


    if(bPatchID() == -1)
    {
        FatalErrorIn("trackedSurface::makeInterpolators()")
            << "Surface patch B not defined."
            << abort(FatalError);
    }

//     patchToPatchInterpolation::setDirectHitTol(1e-2);

    interpolatorBAPtr_ = new IOpatchToPatchInterpolation
    (
        IOobject
        (
            "baInterpolator",
            DB().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh().boundaryMesh()[bPatchID()],
        mesh().boundaryMesh()[aPatchID()],
        intersection::VISIBLE
        // intersection::HALF_RAY
    );


    const scalarField& faceDistBA =
        interpolatorBAPtr_->faceDistanceToIntersection();

    forAll(faceDistBA, faceI)
    {
        if(mag(faceDistBA[faceI] - GREAT) < SMALL)
        {
            FatalErrorIn("trackedSurface::makeInterpolators()")
                << "Error in B-to-A face patchToPatchInterpolation."
                << abort(FatalError);
        }
    }

    const scalarField& pointDistBA =
        interpolatorBAPtr_->pointDistanceToIntersection();

    forAll(pointDistBA, pointI)
    {
        if(mag(pointDistBA[pointI] - GREAT) < SMALL)
        {
            FatalErrorIn("trackedSurface::makeInterpolators()")
                << "Error in B-to-A point patchToPatchInterpolation."
                << abort(FatalError);
        }
    }


    interpolatorABPtr_ = new IOpatchToPatchInterpolation
    (
        IOobject
        (
            "abInterpolator",
            DB().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh().boundaryMesh()[aPatchID()],
        mesh().boundaryMesh()[bPatchID()],
        intersection::VISIBLE
        // intersection::HALF_RAY
    );


    const scalarField& faceDistAB =
        interpolatorABPtr_->faceDistanceToIntersection();

    forAll(faceDistAB, faceI)
    {
        if(mag(faceDistAB[faceI] - GREAT) < SMALL)
        {
            FatalErrorIn("trackedSurface::makeInterpolators()")
                << "Error in A-to-B face patchToPatchInterpolation."
                << abort(FatalError);
        }
    }

    const scalarField& pointDistAB =
        interpolatorABPtr_->pointDistanceToIntersection();

    forAll(pointDistAB, pointI)
    {
        if(mag(pointDistAB[pointI] - GREAT) < SMALL)
        {
            FatalErrorIn("trackedSurface::makeInterpolators()")
                << "Error in A-to-B point patchToPatchInterpolation."
                << abort(FatalError);
        }
    }


    Info << "\nCheck A-to-B and B-to-A interpolators" << endl;

    scalar maxDist = max
    (
        mag
        (
            interpolatorABPtr_->faceInterpolate
            (
                vectorField(mesh().boundaryMesh()[aPatchID()]
               .faceCentres())
            )
          - mesh().boundaryMesh()[bPatchID()].faceCentres()
        )
    );

    scalar maxDistPt = max
    (
        mag
        (
            interpolatorABPtr_->pointInterpolate
            (
                vectorField(mesh().boundaryMesh()[aPatchID()]
               .localPoints())
            )
          - mesh().boundaryMesh()[bPatchID()].localPoints()
        )
    );

    Info << "A-to-B interpolation error, face: " << maxDist
        << ", point: " << maxDistPt << endl;


    maxDist = max
    (
        mag
        (
            interpolatorBAPtr_->faceInterpolate
            (
                vectorField
                (
                    mesh().boundaryMesh()[bPatchID()].faceCentres()
                )
            )
          - mesh().boundaryMesh()[aPatchID()].faceCentres()
        )
    );

    maxDistPt = max
    (
        mag
        (
            interpolatorBAPtr_->pointInterpolate
            (
                vectorField
                (
                    mesh().boundaryMesh()[bPatchID()].localPoints()
                )
            )
          - mesh().boundaryMesh()[aPatchID()].localPoints()
        )
    );

    Info << "B-to-A interpolation error, face: " << maxDist
        << ", point: " << maxDistPt << endl;
}


void trackedSurface::makeControlPoints()
{
    if (debug)
    {
        Info<< "trackedSurface::makeControlPoints() : "
            << "Making control points"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (controlPointsPtr_)
    {
        FatalErrorIn("trackedSurface::makeInterpolators()")
            << "patch to patch interpolators already exists"
            << abort(FatalError);
    }

    IOobject controlPointsHeader
    (
        "controlPoints",
        DB().timeName(),
        mesh(),
        IOobject::MUST_READ
    );

    if (controlPointsHeader.headerOk())
    {
        controlPointsPtr_ =
            new vectorIOField
            (
                IOobject
                (
                    "controlPoints",
                    DB().timeName(),
                    mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            );
    }
    else
    {
        controlPointsPtr_ =
            new vectorIOField
            (
                IOobject
                (
                    "controlPoints",
                    DB().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                aMesh().areaCentres().internalField()
            );

        initializeControlPointsPosition();
    }
}


void trackedSurface::makeMotionPointsMask()
{
    if (debug)
    {
        Info<< "trackedSurface::makeMotionPointsMask() : "
            << "Making motion points mask"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (motionPointsMaskPtr_)
    {
        FatalErrorIn("trackedSurface::motionPointsMask()")
            << "motion points mask already exists"
            << abort(FatalError);
    }


    if(aPatchID() == -1)
    {
        FatalErrorIn("trackedSurface::makeMotionPointsMask()")
            << "Surface patch A not defined."
            << abort(FatalError);
    }


    motionPointsMaskPtr_ = new labelList
    (
        mesh().boundaryMesh()[aPatchID()].nPoints(),
        1
    );
}


void trackedSurface::makeDirections()
{
    if (debug)
    {
        Info<< "trackedSurface::makeDirections() : "
            << "Making displacement directions for points and "
            << "control points"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if
    (
        pointsDisplacementDirPtr_ ||
        facesDisplacementDirPtr_
    )
    {
        FatalErrorIn("trackedSurface::makeDirections()")
            << "points and control points displacement directions "
            << "already exists"
            << abort(FatalError);
    }


    if(aPatchID() == -1)
    {
        FatalErrorIn("trackedSurface::makeDirections()")
            << "Surface patch A not defined."
            << abort(FatalError);
    }


    pointsDisplacementDirPtr_ =
        new vectorField
        (
            mesh().boundaryMesh()[aPatchID()].nPoints(),
            vector::zero
        );

    facesDisplacementDirPtr_ =
        new vectorField
        (
            mesh().boundaryMesh()[aPatchID()].size(),
            vector::zero
        );

    if(!normalMotionDir())
    {
        if(mag(motionDir_) < SMALL)
        {
            FatalErrorIn("trackedSurface::makeDirections()")
                << "Zero motion direction"
                    << abort(FatalError);
        }

        facesDisplacementDir() = motionDir_;
        pointsDisplacementDir() = motionDir_;
    }

    updateDisplacementDirections();
}


void trackedSurface::makeTotalDisplacement()
{
    if (debug)
    {
        Info<< "trackedSurface::makeTotalDisplacement() : "
            << "Making zero total points displacement"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (totalDisplacementPtr_)
    {
        FatalErrorIn("trackedSurface::makeTotalDisplacement()")
            << "total points displacement already exists"
            << abort(FatalError);
    }

    totalDisplacementPtr_ =
        new vectorIOField
        (
            IOobject
            (
                "totalDisplacement",
                DB().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            vectorField
            (
                mesh().boundaryMesh()[aPatchID()].nPoints(),
                vector::zero
            )
        );
}


void trackedSurface::readTotalDisplacement()
{
    if (debug)
    {
        Info<< "trackedSurface::readTotalDisplacement() : "
            << "reading total points displacement if present"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (totalDisplacementPtr_)
    {
        FatalErrorIn("trackedSurface::makeTotalDisplacement()")
            << "total points displacement already exists"
            << abort(FatalError);
    }

    if
    (
        IOobject
        (
            "totalDisplacement",
            DB().timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ).headerOk()
    )
    {
        totalDisplacementPtr_ =
            new vectorIOField
            (
                IOobject
                (
                    "totalDisplacement",
                    DB().timeName(),
                    mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            );
    }
}


void trackedSurface::makeFaMesh() const
{
    if (debug)
    {
        Info<< "trackedSurface::makeFaMesh() : "
            << "Making finite area mesh"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (aMeshPtr_)
    {
        FatalErrorIn("trackedSurface::makeFaMesh()")
            << "finite area mesh already exists"
            << abort(FatalError);
    }

    aMeshPtr_ = new faMesh(mesh());
}

void trackedSurface::makeUs() const
{
    if (debug)
    {
        Info<< "trackedSurface::makeUs() : "
            << "Making surface velocity field"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (UsPtr_)
    {
        FatalErrorIn("trackedSurface::makeUs()")
            << "Surface velocity field already exists"
            << abort(FatalError);
    }


    wordList patchFieldTypes
    (
        aMesh().boundary().size(),
        zeroGradientFaPatchVectorField::typeName
    );

    forAll(aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type()
         == wedgeFaPatch::typeName
        )
        {
            patchFieldTypes[patchI] =
                wedgeFaPatchVectorField::typeName;
        }
        else
        {
            label ngbPolyPatchID =
                aMesh().boundary()[patchI].ngbPolyPatchIndex();

            if (ngbPolyPatchID != -1)
            {
                if
                (
                    mesh().boundary()[ngbPolyPatchID].type()
                 == wallFvPatch::typeName
                )
                {
                    patchFieldTypes[patchI] =
                        slipFaPatchVectorField::typeName;
                }
            }
        }
    }


    UsPtr_ = new areaVectorField
    (
        IOobject
        (
            "Us",
            DB().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aMesh(),
        dimensioned<vector>("Us", dimVelocity, vector::zero),
        patchFieldTypes
    );
}


void trackedSurface::makePhis()
{
    if (debug)
    {
        Info<< "trackedSurface::makePhis() : "
            << "Making surface fluid flux"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (phisPtr_)
    {
        FatalErrorIn("trackedSurface::makePhis()")
            << "Surface fluid flux already exists"
            << abort(FatalError);
    }


    phisPtr_ = new edgeScalarField
    (
        IOobject
        (
            "phis",
            DB().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearEdgeInterpolate(Us()) & aMesh().Le()
    );
}


void trackedSurface::makeSurfactConc() const
{
    if (debug)
    {
        Info<< "trackedSurface::makeSurfactConc() : "
            << "Making surface surfactant concentration field"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (surfactConcPtr_)
    {
        FatalErrorIn("trackedSurface::makeSurfaceConc()")
            << "Surface surfactant concentratio field already exists"
            << abort(FatalError);
    }

    surfactConcPtr_ = new areaScalarField
    (
        IOobject
        (
            "Cs",
            DB().timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh()
    );
}


void trackedSurface::makeSurfaceTension() const
{
    if (debug)
    {
        Info<< "trackedSurface::makeSurfaceTension() : "
            << "Making surface tension field"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (surfaceTensionPtr_)
    {
        FatalErrorIn("trackedSurface::makeSurfaceTension()")
            << "Surface tension field already exists"
            << abort(FatalError);
    }


    surfaceTensionPtr_ = new areaScalarField
    (
        IOobject
        (
            "surfaceTension",
            DB().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cleanInterfaceSurfTension()
      + surfactant().surfactR()*
        surfactant().surfactT()*
        surfactant().surfactSaturatedConc()*
        log(1.0 - surfactantConcentration()/
        surfactant().surfactSaturatedConc())
    );
}


void trackedSurface::makeSurfactant() const
{
    if (debug)
    {
        Info<< "trackedSurface::makeSurfactant() : "
            << "Making surfactant properties"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (surfactantPtr_)
    {
        FatalErrorIn("trackedSurface::makeSurfactant()")
            << "surfactant properties already exists"
            << abort(FatalError);
    }


    const dictionary& surfactProp =
        this->subDict("surfactantProperties");

    surfactantPtr_ = new surfactantProperties(surfactProp);
}


void trackedSurface::makeFluidIndicator()
{
    if (debug)
    {
        Info<< "trackedSurface::makeFluidIndicator() : "
            << "Making fluid indicator"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (fluidIndicatorPtr_)
    {
        FatalErrorIn("trackedSurface::makeFluidIndicator()")
            << "fluid indicator already exists"
            << abort(FatalError);
    }

    fluidIndicatorPtr_ = new volScalarField
    (
        IOobject
        (
            "fluidIndicator",
            DB().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("1", dimless, 1.0),
        zeroGradientFvPatchScalarField::typeName
    );

    volScalarField& fluidIndicator = *fluidIndicatorPtr_;

    if (twoFluids())
    {
        // find start cell
        label pointOnShadowPatch =
            mesh().boundaryMesh()[bPatchID()][0][0];

        label startCell = mesh().pointCells()[pointOnShadowPatch][0];


        // get cell-cells addressing
        const labelListList& cellCells = mesh().cellCells();

        SLList<label> slList(startCell);

        while (slList.size())
        {
            label curCell = slList.removeHead();

            if (fluidIndicator[curCell] == 1)
            {
                fluidIndicator[curCell] = 0.0;

                for (int i = 0; i < cellCells[curCell].size(); i++)
                {
                    slList.append(cellCells[curCell][i]);
                }
            }
        }
    }

    fluidIndicator.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const IOpatchToPatchInterpolation& trackedSurface::interpolatorAB()
{
    if (!interpolatorABPtr_)
    {
        makeInterpolators();
    }

    return *interpolatorABPtr_;
}


const IOpatchToPatchInterpolation& trackedSurface::interpolatorBA()
{
    if (!interpolatorBAPtr_)
    {
        makeInterpolators();
    }

    return *interpolatorBAPtr_;
}


vectorField& trackedSurface::controlPoints()
{
    if (!controlPointsPtr_)
    {
        makeControlPoints();
    }

    return *controlPointsPtr_;
}


labelList& trackedSurface::motionPointsMask()
{
    if (!motionPointsMaskPtr_)
    {
        makeMotionPointsMask();
    }

    return *motionPointsMaskPtr_;
}


vectorField& trackedSurface::pointsDisplacementDir()
{
    if (!pointsDisplacementDirPtr_)
    {
        makeDirections();
    }

    return *pointsDisplacementDirPtr_;
}


vectorField& trackedSurface::facesDisplacementDir()
{
    if (!facesDisplacementDirPtr_)
    {
        makeDirections();
    }

    return *facesDisplacementDirPtr_;
}


vectorField& trackedSurface::totalDisplacement()
{
    if (!totalDisplacementPtr_)
    {
        makeTotalDisplacement();
    }

    return *totalDisplacementPtr_;
}


faMesh& trackedSurface::aMesh()
{
    if (!aMeshPtr_)
    {
        makeFaMesh();
    }

    return *aMeshPtr_;
}

const faMesh& trackedSurface::aMesh() const
{
    if (!aMeshPtr_)
    {
        makeFaMesh();
    }

    return *aMeshPtr_;
}

areaVectorField& trackedSurface::Us()
{
    if (!UsPtr_)
    {
        makeUs();
    }

    return *UsPtr_;
}

const areaVectorField& trackedSurface::Us() const
{
    if (!UsPtr_)
    {
        makeUs();
    }

    return *UsPtr_;
}

edgeScalarField& trackedSurface::Phis()
{
    if (!phisPtr_)
    {
        makePhis();
    }

    return *phisPtr_;
}

areaScalarField& trackedSurface::surfactantConcentration()
{
    if (!surfactConcPtr_)
    {
        makeSurfactConc();
    }

    return *surfactConcPtr_;
}

const areaScalarField& trackedSurface::surfactantConcentration() const
{
    if (!surfactConcPtr_)
    {
        makeSurfactConc();
    }

    return *surfactConcPtr_;
}

areaScalarField& trackedSurface::surfaceTension()
{
    if (!surfaceTensionPtr_)
    {
        makeSurfaceTension();
    }

    return *surfaceTensionPtr_;
}

const areaScalarField& trackedSurface::surfaceTension() const
{
    if (!surfaceTensionPtr_)
    {
        makeSurfaceTension();
    }

    return *surfaceTensionPtr_;
}

const surfactantProperties& trackedSurface::surfactant() const
{
    if (!surfactantPtr_)
    {
        makeSurfactant();
    }

    return *surfactantPtr_;
}


const volScalarField& trackedSurface::fluidIndicator()
{
    if (!fluidIndicatorPtr_)
    {
        makeFluidIndicator();
    }

    return *fluidIndicatorPtr_;
}


tmp<areaVectorField> trackedSurface::surfaceTensionGrad()
{
    tmp<areaVectorField> tgrad
    (
        new areaVectorField
        (
            IOobject
            (
                "surfaceTensionGrad",
                DB().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            (-fac::grad(surfactantConcentration())*
            surfactant().surfactR()*surfactant().surfactT()/
            (1.0 - surfactantConcentration()/
            surfactant().surfactSaturatedConc()))()
        )
    );

    return tgrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
