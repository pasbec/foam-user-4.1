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
#include "emptyFaPatch.H"
#include "wedgeFaPatch.H"
#include "PstreamCombineReduceOps.H"
#include "coordinateSystem.H"
#include "scalarMatrices.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void trackedSurface::writeVTK() const
{
    aMesh().patch().writeVTK
    (
        DB().timePath()/"trackedSurface",
        aMesh().patch(),
        aMesh().patch().points()
    );
}


void trackedSurface::writeVTKControlPoints()
{
    // Write patch and points into VTK
    fileName name(DB().timePath()/"trackedSurfaceControlPoints");
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
    const label aPatchID,
    const fvMesh& vMesh,
    const GeometricField<Type, faPatchField, areaMesh>& af
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "vol" + af.name(),
                af.instance(),
                af.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            vMesh,
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

    forAll(vMesh.boundaryMesh(), patchI)
    {
        vf.boundaryField()[patchI] == pTraits<Type>::zero;
    }

    vf.boundaryField()[aPatchID] == af;

    vf.write();

    tvf.clear();
}


void trackedSurface::writeVolA()
{
    Us().write();

    (fac::div(Us()))().write();

    (fac::grad(Us()))().write();

    aMesh().Le().write();

    (fac::average(aMesh().Le()))().write();

    aMesh().edgeAreaNormals().write();

    aMesh().faceAreaNormals().write();

    aMesh().faceCurvatures().write();

    surfaceTensionGrad()().write();

    // Us
    writeVol(aPatchID(), mesh(), Us());

    // fac::div(Us())
    writeVol(aPatchID(), mesh(), (fac::div(Us()))());

    // fac::grad(Us())
    writeVol(aPatchID(), mesh(), (fac::grad(Us()))());

    // volLeAverage
    writeVol(aPatchID(), mesh(), (fac::average(aMesh().Le()))());

    // faceAreaNormals
    writeVol(aPatchID(), mesh(), aMesh().faceAreaNormals());

    // faceCurvatures
    writeVol(aPatchID(), mesh(), aMesh().faceCurvatures());

    // faceCurvaturesDivNormals
    writeVol
    (
        aPatchID(),
        mesh(),
        areaScalarField
        (
            "faceCurvaturesDivNormals",
            -fac::div(aMesh().faceAreaNormals())
        )
    );

    // surfaceTensionGrad
    writeVol(aPatchID(), mesh(), surfaceTensionGrad()());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
