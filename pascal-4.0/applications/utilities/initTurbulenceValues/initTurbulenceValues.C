#include "fvCFD.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    // Add region option
    #include "addRegionOption.H"

    // Add time selector options
    timeSelector::addOptions();

    // Add export file name option
    argList::validOptions.insert("file", "name");

    // Set root case and create time
    #include "setRootCase.H"
    #include "createTime.H"

    // Create time instant list
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

    // Info about what is going on
    Info << "Initialize uniform turbulence internal field values" << endl << nl;

    // Loop over all time instances
    forAll(times, timeI)
    {
        // Set time and update mesh
        runTime.setTime(times[timeI], timeI);
        mesh.readUpdate();

        Info << "Time = " << runTime.timeName() << endl << nl;

        // Read turbulence Properties dictionary
        IOdictionary turbulenceProperties
        (
            IOobject
            (
                "turbulenceProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            )
        );

        // Read simulationType
        word simulationType = pTraits<word>(turbulenceProperties.lookup("simulationType"));

        // Read initialization sub dictionary
        dictionary& turbulencePropertiesInitSubDict = turbulenceProperties.subDict("init");

        scalar refU = pTraits<scalar>(turbulencePropertiesInitSubDict.lookup("refU"));
        scalar refI = pTraits<scalar>(turbulencePropertiesInitSubDict.lookup("refI"));
        scalar refL = pTraits<scalar>(turbulencePropertiesInitSubDict.lookup("refL"));

        // List of all "turbulence" fields
        word foList[] = {"k", "epsilon", "omega"};
        int foLostSize = sizeof(foList) / sizeof(foList[0]);

        for (int i=0; i<foLostSize; i++)
        {
            word vfName = foList[i];

            // Read base object if present
            IOobject fo
            (
                vfName,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT
            );

            // Check if field exists
            if (fo.headerOk())
            {
                // Read field
                volScalarField vf
                (
                    IOobject
                    (
                        vfName,
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE,
                        true
                    ),
                    mesh
                );

                // Calculate estimations of turbulence properties based on dimensional analysis
                // http://www.cfd-online.com/Wiki/Turbulence_modeling
                scalar Cmu       = 0.09;
                scalar half      = 0.5;
                scalar threeHalf = 3.0/2.0;

                scalar refK = threeHalf * std::pow(refU*refI, 2);
                scalar refE = Cmu * std::pow(refK, threeHalf) / refL;
                scalar refO = std::pow(refK, half) / refL;

                if (vfName == "k")
                {
                    vf.internalField() = refK;
                    Info << " " << vfName << " = " << refK << endl;
                }
                else if (vfName == "epsilon")
                {
                    vf.internalField() = refE;
                    Info << " " << vfName << " = " << refE << endl;
                }
                else if (vfName == "omega")
                {
                    vf.internalField() = refO;
                    Info << " " << vfName << " = " << refO << endl;
                }

                // Write new values to field
                vf.write();
            }
        }
        runTime.write();
    }

    Info << "\nEnd\n" << endl;

    return(0);
}

