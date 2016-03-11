/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

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

#include "coordinateSystem.H"
#include "scalarMatrices.H"
#include "zeroGradientFaPatchFields.H"
#include "fixedGradientFaPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void trackedSurface::writeVTK() const
{
    aMesh().patch().writeVTK
    (
        DB().timePath()/prefix_,
        aMesh().patch(),
        aMesh().patch().points()
    );
}


void trackedSurface::writeVTKControlPoints()
{
    // Write patch and points into VTK
    fileName name(DB().timePath()/prefix_+"ControlPoints");
    OFstream mps(name + ".vtk");

    mps << "# vtk DataFile Version 2.0" << nl
        << name << ".vtk" << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl
        << "POINTS " << controlPoints().size() << " float" << nl;

    forAll(controlPoints(), pointI)
    {
        mps << controlPoints()[pointI].x() << ' '
            << controlPoints()[pointI].y() << ' '
            << controlPoints()[pointI].z() << nl;
    }

    // Write vertices
    mps << "VERTICES " << controlPoints().size() << ' '
        << controlPoints().size()*2 << nl;

    forAll(controlPoints(), pointI)
    {
        mps << 1 << ' ' << pointI << nl;
    }
}


template<class Type>
void
trackedSurface::writeVol
(
    const GeometricField<Type, faPatchField, areaMesh>& af
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "vol" + word(toupper(af.name()[0]))
                    + word(af.name().substr(1)),
                DB().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh(),
            dimensioned<Type>
            (
                "0",
                af.dimensions(),
                pTraits<Type>::zero
            ),
            fixedValueFvPatchField<Type>::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf();

    forAll(mesh().boundaryMesh(), patchI)
    {
        vf.boundaryField()[patchI] == pTraits<Type>::zero;
    }

    vf.boundaryField()[aPatchID()] == af;

    vf.write();

    tvf.clear();
}


void trackedSurface::writeVolA()
{
    // Us
    writeVol(Us());

    // fac::div(Us())
    writeVol((fac::div(Us()))());

    // fac::grad(Us())
    writeVol((fac::grad(Us()))());

    // faceAreaNormals
    writeVol(aMesh().faceAreaNormals());

    // faceCurvatures
    writeVol(aMesh().faceCurvatures());

    // faceCurvaturesDivNormals
    writeVol
    (
        areaScalarField
        (
            "faceCurvaturesDivNormals",
            -fac::div(aMesh().faceAreaNormals())
        )
    );

    // surfaceTensionGrad
    writeVol(surfaceTensionGrad()());

    // H
    tmp<areaVectorField> tH
    (
        new areaVectorField
        (
            IOobject
            (
                "H",
                DB().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            aMesh(),
            dimensionedVector
            (
                word(),
                dimless,
                vector::zero
            ),
            zeroGradientFaPatchVectorField::typeName
        )
    );
    areaVectorField& H = tH();

    H.internalField() =
        controlPoints() - aMesh().areaCentres().internalField();
    H.correctBoundaryConditions();

    writeVol(H);

    // H
    tmp<areaScalarField> tmagH
    (
        new areaScalarField
        (
            IOobject
            (
                "magH",
                DB().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            aMesh(),
            dimensionedScalar
            (
                word(),
                dimless,
                0
            ),
            zeroGradientFaPatchScalarField::typeName
        )
    );
    areaScalarField& magH = tmagH();

    magH.internalField() =
        mag(H.internalField())
      * sign(H.internalField()&facesDisplacementDir());
    magH.correctBoundaryConditions();

    writeVol(magH);
    tmagH.clear();
    tH.clear();

    // nGradUn
    tmp<areaScalarField> tnGradUn
    (
        new areaScalarField
        (
            IOobject
            (
                "nGradUn",
                DB().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            aMesh(),
            dimensionedScalar
            (
                word(),
                dimless,
                0
            ),
            zeroGradientFaPatchScalarField::typeName
        )
    );
    areaScalarField& anGradUn = tnGradUn();

    anGradUn.internalField() = nGradUn();
    anGradUn.correctBoundaryConditions();

    writeVol(anGradUn);
    tnGradUn.clear();

    // phi
    tmp<areaScalarField> taPhi
    (
        new areaScalarField
        (
            IOobject
            (
                "aPhi",
                DB().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            aMesh(),
            dimensionedScalar
            (
                word(),
                dimless,
                0
            ),
            zeroGradientFaPatchScalarField::typeName
        )
    );
    areaScalarField& aPhi = taPhi();

    aPhi.internalField() =
        phi_.boundaryField()[aPatchID()];
    aPhi.correctBoundaryConditions();

    writeVol(aPhi);

    // meshPhi
    tmp<areaScalarField> taPhiMesh
    (
        new areaScalarField
        (
            IOobject
            (
                "aPhiMesh",
                DB().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            aMesh(),
            dimensionedScalar
            (
                word(),
                dimless,
                0
            ),
            zeroGradientFaPatchScalarField::typeName
        )
    );
    areaScalarField& aPhiMesh = taPhiMesh();

    aPhiMesh.internalField() =
        fvc::meshPhi(rho(),U())().boundaryField()[aPatchID()];
    aPhiMesh.correctBoundaryConditions();

    writeVol(aPhiMesh);

    // diffPhi
    writeVol
    (
        areaScalarField
        (
            "aPhiDiff",
            aPhi-aPhiMesh
        )
    );
    taPhi.clear();
    taPhiMesh.clear();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
