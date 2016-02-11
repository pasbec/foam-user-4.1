
#include "argList.H"
#include "foamTime.H"
#include "polyMesh.H"
#include "meshTools.H"

#include "cellZone.H"
#include "faceZone.H"
#include "pointZone.H"

#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"


using namespace Foam;

int main(int argc, char *argv[])
{
#   include "addRegionOption.H"

#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();

    Foam::word meshRegionName = polyMesh::defaultRegion;
    args.optionReadIfPresent("region", meshRegionName);

#   include "createNamedPolyMesh.H"



    // cellZones

    Info << "cellZones:" << endl
        << mesh.cellZones().size() << endl
        << "(" << endl;

    forAll(mesh.cellZones(), cellZoneI)
    {
        cellZone& curCellZone = mesh.cellZones()[cellZoneI];

        Info << "    " << curCellZone.name() << endl;

        labelHashSet curCellZoneLabelHashSet(curCellZone);

        cellSet curCellSet(mesh, "cellSet_"+curCellZone.name(), curCellZoneLabelHashSet);

        curCellSet.Foam::regIOobject::write();
    }

    Info << ")" << endl << nl;



    // faceZones

    Info << "faceZones:" << endl
        << mesh.faceZones().size() << endl
        << "(" << endl;

    forAll(mesh.faceZones(), faceZoneI)
    {
        faceZone& curFaceZone = mesh.faceZones()[faceZoneI];

        Info << "    " << curFaceZone.name() << endl;

        labelHashSet curFaceZoneLabelHashSet(curFaceZone);

        faceSet curFaceSet(mesh, "faceSet_"+curFaceZone.name(), curFaceZoneLabelHashSet);

        curFaceSet.Foam::regIOobject::write();
    }

    Info << ")" << endl << nl;



    // pointZones

    Info << "pointZones:" << endl
        << mesh.pointZones().size() << endl
        << "(" << endl;

    forAll(mesh.pointZones(), pointZoneI)
    {
        pointZone& curPointZone = mesh.pointZones()[pointZoneI];

        Info << "    " << curPointZone.name() << endl;

        labelHashSet curPointZoneLabelHashSet(curPointZone);

        pointSet curPointSet(mesh, "pointSet_"+curPointZone.name(), curPointZoneLabelHashSet);

        curPointSet.Foam::regIOobject::write();
    }

    Info << ")" << endl << nl;



    Info << "\nEnd\n" << endl;

    return(0);
}

