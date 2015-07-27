#include "fvCFD.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    // Add region option
    #include "addRegionOption.H"

    // Add command line options for time instant selection
    timeSelector::addOptions();

    // Add command line argument to provide a patch name
    argList::validArgs.append("patchName"); // arg #1

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

        // Try to find patch label for given patch name from command line arg #1
        const word exPatchName (IStringStream(args.args()[1])());
        const label exPatchID = mesh.boundaryMesh().findPatchID(exPatchName);

        if (exPatchID > -1)
        {

            // Reference to the extraction data
            const vectorField& exData = mesh.boundaryMesh()[exPatchID].localPoints();

            // Set file names
            fileName exDir(runTime.timeName()/"rawdata");
            fileName exFile(exDir/"dataPoints."+exPatchName+".dat");

            // Print info
            Info << "Extract patch point data to file: " << exFile << endl << endl;

            // Create output directory and stream
            mkDir(exDir);
            OFstream ex(exFile);

            ex << "# Source = Patch, ID: " << exPatchID << ", Name: " << exPatchName << endl;
            ex << "# Time   = " << runTime.timeName() << endl;
            ex << "# Data   = Mesh point coordinates (1:x 2:y 3:z)" << endl;

            // Extract all data
            forAll(exData, exDataI)
            {

                ex << exData[exDataI].component(0) << " "
                   << exData[exDataI].component(1) << " "
                   << exData[exDataI].component(2) << endl;

            }

        } else {

              FatalError
                  << "A patch called '" << exPatchName << "' does not exist!"
                  << abort(FatalError);

        }
    }

    Info << "End" << endl << endl;

}

