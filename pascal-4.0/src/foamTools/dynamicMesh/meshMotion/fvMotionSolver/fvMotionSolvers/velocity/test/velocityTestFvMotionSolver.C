/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

#include "velocityTestFvMotionSolver.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "volPointInterpolation.H"
#include "leastSquaresVolPointInterpolation.H"
#include "fvc.H"
#include "slipFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(velocityTestFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        fvMotionSolver,
        velocityTestFvMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::velocityTestFvMotionSolver::velocityTestFvMotionSolver
(
    const polyMesh& mesh,
    Istream&
)
:
    fvMotionSolver(mesh),
    pointMotionU_
    (
        IOobject
        (
            "pointMotionU",
            fvMesh_.time().timeName(),
            fvMesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pointMesh::New(fvMesh_)
    ),
    cellMotionU_
    (
        IOobject
        (
            "cellMotionU",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector
        (
            "cellMotionU",
            pointMotionU_.dimensions(),
            vector::zero
        ),
        cellMotionBoundaryTypes<vector>(pointMotionU_.boundaryField())
    ),
    diffusivityPtr_
    (
        motionDiffusivity::New(*this, lookup("diffusivity"))
    ),
    leastSquaresVolPoint_(false)
{
    if (found("leastSquaresVolPoint"))
    {
        leastSquaresVolPoint_ =
            Switch(lookup("leastSquaresVolPoint"));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::velocityTestFvMotionSolver::~velocityTestFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::velocityTestFvMotionSolver::curPoints() const
{
    if (leastSquaresVolPoint_)
    {
        leastSquaresVolPointInterpolation::New(fvMesh_).interpolate
        (
            cellMotionU_,
            pointMotionU_
        );
    }
    else
    {
        volPointInterpolation::New(fvMesh_).interpolate
        (
            cellMotionU_,
            pointMotionU_
        );
    }

    tmp<pointField> tcurPoints(new pointField(fvMesh_.allPoints()));
    pointField& cp = tcurPoints();
    const pointField& pointMotionUI = pointMotionU_.internalField();

    forAll(pointMotionUI, pointI)
    {
        cp[pointI] +=
            pointMotionUI[pointI]*fvMesh_.time().deltaT().value();
    }

    // tmp<pointField> tcurPoints
    // (
    //     fvMesh_.points()
    //   + fvMesh_.time().deltaT().value()*pointMotionU_.internalField()
    // );

    twoDCorrectPoints(tcurPoints());

    return tcurPoints;
}


void Foam::velocityTestFvMotionSolver::solve()
{
    // The points have moved so before interpolation update
    // the fvMotionSolver accordingly
    movePoints(fvMesh_.allPoints());

    diffusivityPtr_->correct();

    // ZT, Problem on symmetry plane
//     pointMotionU_.boundaryField().updateCoeffs();

// TEST

    volScalarField motionDiffusivity("motionDiffusivity", fvc::average(diffusivityPtr_->operator()()));
    motionDiffusivity.write();

// TEST

    surfaceScalarField meshDelta
    (
        IOobject
        (
            "meshDelta",
            fvMesh_.time().timeName(),
            fvMesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        1.0/fvMesh_.deltaCoeffs()
    );

    IOobject refMeshDeltaObj
    (
        "meshDelta",
        fvMesh_.time().constant(),
        fvMesh_,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    surfaceScalarField refMeshDelta
    (
        refMeshDeltaObj,
        meshDelta
    );

    if (!refMeshDeltaObj.headerOk())
    {
        refMeshDelta.write();
    }

// TEST

    label nNonOrthCorr = 1;
    if (found("nNonOrthogonalCorrectors"))
    {
        nNonOrthCorr = readInt(lookup("nNonOrthogonalCorrectors"));
    }

    for (label i=0; i<nNonOrthCorr; i++)
    {
        Pout << "Correction: " << i << endl;

        scalar nu_ = 0.3;

        surfaceScalarField& muf = diffusivityPtr_->operator()()();
        surfaceScalarField lambdaf(word(), muf*(2*nu_/(1 - 2*nu_)));

        volScalarField mu
        (
            "mu",
            fvc::average(muf)
          / dimensionedScalar(word(), muf.dimensions(), 1.0)
        );

        volScalarField lambda
        (
            "lambda",
            fvc::average(lambdaf)
          / dimensionedScalar(word(), lambdaf.dimensions(), 1.0)
        );

        volVectorField& U = cellMotionU_;
        volTensorField gradU("gradCellMotionU", fvc::grad(U));

        wordList FpatchTypes(fvMesh_.boundary().size());

        forAll (U.boundaryField(), patchI)
        {
            if (isA<slipFvPatchVectorField>(U.boundaryField()[patchI]))
            {
                FpatchTypes[patchI] = zeroGradientFvPatchVectorField::typeName;
            }
            else if (isA<cellMotionFvPatchVectorField>(U.boundaryField()[patchI]))
            {
                FpatchTypes[patchI] = fixedValueFvPatchVectorField::typeName;
            }
            else
            {
                FpatchTypes[patchI] = U.boundaryField()[patchI].type();
            }
        }

        volVectorField F
        (
            IOobject
            (
                "cellMotionF",
                fvMesh_.time().timeName(),
                fvMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvMesh_,
            dimensionedVector(word(), dimless/dimLength/dimTime, vector::zero),
            FpatchTypes
        );

        forAll (F.boundaryField(), patchI)
        {
            fvPatchVectorField& pf = F.boundaryField()[patchI];

            if (isA<fixedValueFvPatchVectorField>(pf))
            {
                F.boundaryField()[patchI] ==
                    1.0
                  * (meshDelta.boundaryField()[patchI] - refMeshDelta.boundaryField()[patchI]) / refMeshDelta.boundaryField()[patchI]
                  / refMeshDelta.boundaryField()[patchI]
                  / fvMesh_.time().deltaT().value()
                  * fvMesh_.boundary()[patchI].nf();
            }
        }

        Foam::solve
        (
            fvm::laplacian
            (
                F
            )
        );

//         F = U / (mag(U)+dimensionedScalar(word(), U.dimensions(), VSMALL)) / dimensionedScalar(word(), dimLength*dimTime, 1.0);

//         F = (F/dimensionedScalar(word(), F.dimensions(), 1.0) & gradRelDiffMeshV/dimensionedScalar(word(), gradRelDiffMeshV.dimensions(), 1.0)) * F;

//         F.replace(0, 1e+6*relDiffMeshDeltaX);
//         F.replace(1, 1e+6*relDiffMeshDeltaY);
//         F.replace(2, 1e+6*relDiffMeshDeltaZ);

//         F.replace(0, gradXRelDiffMeshDeltaX);
//         F.replace(1, gradYRelDiffMeshDeltaY);
//         F.replace(2, gradZRelDiffMeshDeltaZ);

        F *= 1e-7;
        F.internalField() /= fvMesh_.V();
        F.correctBoundaryConditions();

        F.write();

        Foam::solve
        (
            fvm::laplacian
            (
                2*mu + lambda, U,
                "laplacian(diffusivity,cellMotionU)"
            )
          + fvc::div
            (
                mu*gradU.T() + lambda*(I*tr(gradU)) - (mu + lambda)*gradU,
                "div(cellMotionSigma)"
            )
          + F
        );
    }

// TEST

//     volVectorField smoothCellMotionU
//     (
//         IOobject
//         (
//             "cellMotionU",
//             fvMesh_.time().timeName(),
//             fvMesh_,
//             IOobject::NO_READ,
//             IOobject::AUTO_WRITE,
//             false
//         ),
//         fvMesh_,
//         dimensionedVector(word(), cellMotionU_.dimensions(), vector::zero),
//         slipFvPatchVectorField::typeName
//     );
//
//     volVectorField smoothSource
//     (
//         IOobject
//         (
//             "smoothSource",
//             fvMesh_.time().timeName(),
//             fvMesh_,
//             IOobject::NO_READ,
//             IOobject::AUTO_WRITE,
//             false
//         ),
//         fvMesh_,
//         dimensionedVector(word(), cellMotionU_.dimensions()/dimArea, vector::zero),
//         zeroGradientFvPatchVectorField::typeName
//     );
//
// //     smoothSource.replace(0, -relDiffMeshDeltaX);
// //     smoothSource.replace(1, -relDiffMeshDeltaY);
// //     smoothSource.replace(2, -relDiffMeshDeltaZ);
//
// //     smoothSource.replace(0, -(relDiffMeshV-gAverage(relDiffMeshV)));
// //     smoothSource.replace(1, -(relDiffMeshV-gAverage(relDiffMeshV)));
// //     smoothSource.replace(2, -(relDiffMeshV-gAverage(relDiffMeshV)));
//
// //     smoothSource.replace(0, gradRelDiffMeshV.component(0)-average(gradRelDiffMeshV.component(0)));
// //     smoothSource.replace(1, gradRelDiffMeshV.component(1)-average(gradRelDiffMeshV.component(1)));
// //     smoothSource.replace(2, gradRelDiffMeshV.component(2)-average(gradRelDiffMeshV.component(2)));
//
//     smoothSource.write();
//
//     Foam::solve
//     (
//         fvm::laplacian
//         (
//             smoothCellMotionU,
//             "laplacian(diffusivity,cellMotionU)"
//         )
//      == smoothSource
//     );
//
//     smoothCellMotionU.rename("smoothCellMotionU");
//     smoothCellMotionU.write();

// TEST

//     cellMotionU_ == cellMotionU_ + smoothCellMotionU;

}


void Foam::velocityTestFvMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    fvMotionSolver::updateMesh(mpm);

    // Update diffusivity. Note two stage to make sure old one is de-registered
    // before creating/registering new one.
    diffusivityPtr_.reset(NULL);
    diffusivityPtr_ = motionDiffusivity::New(*this, lookup("diffusivity"));
}


// ************************************************************************* //
