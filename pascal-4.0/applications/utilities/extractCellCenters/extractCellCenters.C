#include "fvCFD.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    // Add region option
    #include "addRegionOption.H"

    // Add command line options for time instant selection
    timeSelector::addOptions();

    // Add export file name option
    argList::validOptions.insert("file", "name");

#   include "setRootCase.H"
#   include "createTime.H"

    instantList times = timeSelector::select0(runTime, args);

    // Read region name
    Foam::word meshRegionName = polyMesh::defaultRegion;
    args.optionReadIfPresent("region", meshRegionName);

    // Create Mesh
    fvMesh mesh
    (
        IOobject
        (
            meshRegionName,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    // Loop over all times from instant list
    forAll(times, timeI)
    {

        // Set time and update mesh
        runTime.setTime(times[timeI], timeI);
        mesh.readUpdate();

        // Set file names
        fileName exDir(runTime.timeName()/"rawdata");
        //fileName exFile(exDir/argv[0]+".dat");
        fileName exFile(exDir/"dataCenters.internal.dat");

        // Print info
        Info << "Extract cell center data to file: " << exFile << endl << endl;

        // Create output directory and stream
        mkDir(exDir);
        OFstream ex(exFile);

        ex << "# Source = Domain" << endl;
        ex << "# Time   = " << runTime.timeName() << endl;
        ex << "# Data   = Cell center coordinates (1:x 2:y 3:z)" << endl;

        // Reference to the extraction data
        const vectorField& exData = mesh.C();

        // Extract all data
        forAll(exData, exDataI)
        {

            ex << exData[exDataI].component(0) << " "
               << exData[exDataI].component(1) << " "
               << exData[exDataI].component(2) << endl;

        }
    }

    Info << "End" << endl << endl;

}

