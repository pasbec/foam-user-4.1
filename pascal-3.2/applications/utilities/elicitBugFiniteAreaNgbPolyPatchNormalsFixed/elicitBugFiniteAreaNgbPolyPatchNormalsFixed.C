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

Application
    elicitBugFiniteAreaNgbPolyPatchNormals

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "faCFD.H"

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    using namespace Foam;
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    faMesh* aMeshPtr = new faMesh(mesh);
    faMesh& aMesh = *aMeshPtr;

    forAll(aMesh.boundary(), patchI)
    {
        label ngbPolyPatchID = aMesh.boundary()[patchI].ngbPolyPatchIndex();

        string ngbPolyPatchName = mesh.boundary()[ngbPolyPatchID].name();
        string ngbPolyPatchType = mesh.boundary()[ngbPolyPatchID].type();

        vectorField ngbPolyPatchNf = mesh.boundary()[ngbPolyPatchID].nf();

        vectorField ngbPolyPatchFaceNormals = aMesh.boundary()[patchI].ngbPolyPatchFaceNormals();
        vectorField ngbPolyPatchPointNormals = aMesh.boundary()[patchI].ngbPolyPatchPointNormals();

        Info << "Information for neighbour polyPatch of finite area boundary edge patch " << aMesh.boundary()[patchI].name() << " with ID " << patchI << ": " << endl;
        Info << "- " << "ID < aMesh.boundary()[patchI].ngbPolyPatchIndex() >: " << ngbPolyPatchID << endl;
        Info << "- " << "Name < mesh.boundary()[ID].name() >: " << ngbPolyPatchName << endl;
        Info << "- " << "Type < mesh.boundary()[ID].type() >: " << ngbPolyPatchType << endl;
        Info << "- " << "Normals < mesh.boundary()[ID].nf() >: " << ngbPolyPatchNf << endl;
        Info << "- " << "PointNormals < aMesh.boundary()[ID].ngbPolyPatchPointNormals() >: " << ngbPolyPatchPointNormals << endl;
        Info << "- " << "FaceNormals < aMesh.boundary()[ID].ngbPolyPatchFaceNormals() >: " << ngbPolyPatchFaceNormals << endl;
        Info << endl;
    }

    Info << "End\n" << endl;

    return(0);
}

// ************************************************************************* //

