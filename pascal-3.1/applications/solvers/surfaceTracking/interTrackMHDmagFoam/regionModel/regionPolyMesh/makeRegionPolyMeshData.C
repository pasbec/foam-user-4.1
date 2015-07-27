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

#include "regionPolyMesh.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void regionPolyMesh::makeCellRegionMap(const label& regionI) const
{
    if (debug)
    {
        Info << "regionPolyMesh::makeCellRegionMap(...) : "
            << "Create region cell map for region "
            << regionName(regionI)
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (cellRegionMap_[regionI])
    {
        FatalErrorIn("regionPolyMesh::makeCellRegionMap(...)")
            << "Region cell map for region "
            << regionName(regionI)
            << " already exists"
            << abort(FatalError);
    }

    // Create region cell map
    cellRegionMap_[regionI] = newRegionMap(regionI,"cell");
}

void regionPolyMesh::makePointRegionMap(const label& regionI) const
{
    if (debug)
    {
        Info << "regionPolyMesh::makePointRegionMap(...) : "
            << "Create region point map for region "
            << regionName(regionI)
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (pointRegionMap_[regionI])
    {
        FatalErrorIn("regionPolyMesh::makePointRegionMap(...)")
            << "Region point map for region "
            << regionName(regionI)
            << " already exists"
            << abort(FatalError);
    }

    // Create region point map
    pointRegionMap_[regionI] = newRegionMap(regionI,"point");
}

void regionPolyMesh::makeCellProcMap(const label& regionI) const
{
    if (debug)
    {
        Info << "regionPolyMesh::makeCellProcMap(...) : "
            << "Create processor cell map for region "
            << regionName(regionI)
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (cellProcMap_[regionI])
    {
        FatalErrorIn("regionPolyMesh::makeCellProcMap(...)")
            << "Processor cell map for region "
            << regionName(regionI)
            << " already exists"
            << abort(FatalError);
    }

    // Create processor cell map
    cellProcMap_[regionI] = newProcMap(regionI,"cell");
}

void regionPolyMesh::makePointProcMap(const label& regionI) const
{
    if (debug)
    {
        Info << "regionPolyMesh::makePointProcMap(...) : "
            << "Create processor point map for region "
            << regionName(regionI)
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (pointProcMap_[regionI])
    {
        FatalErrorIn("regionPolyMesh::makePointProcMap(...)")
            << "Processor point map for region "
            << regionName(regionI)
            << " already exists"
            << abort(FatalError);
    }

    // Create processor point map
    pointProcMap_[regionI] = newProcMap(regionI,"point");
}

void regionPolyMesh::makeCellRegionProcMap(const label& regionI) const
{
    if (debug)
    {
        Info << "regionPolyMesh::makeCellRegionProcMap(...) : "
            << "Create region-processor cell map for region "
            << regionName(regionI)
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (cellRegionProcMap_[regionI])
    {
        FatalErrorIn("regionPolyMesh::makeCellRegionProcMap(...)")
            << "Region-processor cell map for region "
            << regionName(regionI)
            << " already exists"
            << abort(FatalError);
    }

    // Create processor cell map
    cellRegionProcMap_[regionI] = newCellRegionProcMap(regionI);
}

void regionPolyMesh::makePointRegionProcMap(const label& regionI) const
{
    if (debug)
    {
        Info << "regionPolyMesh::makePointRegionProcMap(...) : "
            << "Create region-processor point map for region "
            << regionName(regionI)
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (pointRegionProcMap_[regionI])
    {
        FatalErrorIn("regionPolyMesh::makePointRegionProcMap(...)")
            << "Region-processor point map for region "
            << regionName(regionI)
            << " already exists"
            << abort(FatalError);
    }

    // Create processor point map
    pointRegionProcMap_[regionI] = newPointRegionProcMap(regionI);
}

void regionPolyMesh::makeCellMap(const label& regionI) const
{
    if (debug)
    {
        Info << "regionPolyMesh::makeCellMap(...) : "
            << "Create cell map for region "
            << regionName(regionI)
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (cellMap_[regionI])
    {
        FatalErrorIn("regionPolyMesh::makeCellMap(...)")
            << "Cell map for region "
            << regionName(regionI)
            << " already exists"
            << abort(FatalError);
    }

    // Create cell map
    if (!parallel())
    {
        // Serial
        if (!cellRegionMap_[regionI])
        {
            makeCellRegionMap(regionI);
        }

        cellMap_[regionI] = cellRegionMap_[regionI];
    }
    else
    {
        if (!parallelSplitRegions())
        {
            // Parallel split region-processor
            if (!cellRegionProcMap_[regionI])
            {
                makeCellRegionProcMap(regionI);
            }

            cellMap_[regionI] = cellRegionProcMap_[regionI];
        }
    }
}

void regionPolyMesh::makePointMap(const label& regionI) const
{
    if (debug)
    {
        Info << "regionPolyMesh::makePointMap(...) : "
            << "Create point map for region "
            << regionName(regionI)
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (pointMap_[regionI])
    {
        FatalErrorIn("regionPolyMesh::makePointMap(...)")
            << "Point map for region "
            << regionName(regionI)
            << " already exists"
            << abort(FatalError);
    }

    // Create point map
    if (!parallel())
    {
        // Serial
        if (!pointRegionMap_[regionI])
        {
            makePointRegionMap(regionI);
        }

        pointMap_[regionI] = pointRegionMap_[regionI];
    }
    else
    {
        if (!parallelSplitRegions())
        {
            // Parallel split region-processor
            if (!pointRegionProcMap_[regionI])
            {
                makePointRegionProcMap(regionI);
            }

            pointMap_[regionI] = pointRegionProcMap_[regionI];
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

labelIOList* regionPolyMesh::newRegionMap
(
    const label& regionI,
    const word& topoType
) const
{
    fileName mapInstance;

    if (!parallel())
    {
        // Serial
        mapInstance = mesh(regionI).time().constant();
    }
    else
    {
        if (!parallelSplitRegions())
        {
            // Parallel split region-processor
            mapInstance = mesh(regionI).time().caseConstant();
        }
    }

    IOobject mapObj
    (
        topoType + word("RegionAddressing"),
        mapInstance,
        polyMesh::meshSubDir,
        mesh(regionI),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    if(regionName(regionI) == polyMesh::defaultRegion)
    {
        return new labelIOList
        (
            mapObj,
            labelList()
        );
    }
    else
    {
        mapObj.readOpt() = IOobject::MUST_READ;
        mapObj.writeOpt() = IOobject::AUTO_WRITE;

        return new labelIOList
        (
            mapObj
        );
    }
}

labelIOList* regionPolyMesh::newProcMap
(
    const label& regionI,
    const word& topoType
) const
{
    word searchName = (topoType == "point") ? "points" : "faces";

    fileName mapInstance;

    if (!parallel())
    {
        // Serial
        mapInstance = mesh(regionI).time().constant();
    }
    else
    {
        if (!parallelSplitRegions())
        {
            // Parallel split region-processor
            mapInstance = mesh(regionI).time().constant();
        }
    }

    IOobject mapObj
    (
        topoType + word("ProcAddressing"),
        mapInstance,
        polyMesh::meshSubDir,
        mesh(regionI),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    return new labelIOList
    (
        mapObj
    );
}

labelIOList* regionPolyMesh::newCellRegionProcMap
(
    const label& regionI
) const
{

    IOobject mapObj
    (
        "cellRegionProcAddressing",
        mesh(regionI).time().constant(),
        polyMesh::meshSubDir,
        mesh(regionI),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    if (regionName(regionI) == polyMesh::defaultRegion)
    {
        return new labelIOList
        (
            mapObj,
            labelList()
        );
    }
    else
    {
        mapObj.readOpt() = IOobject::MUST_READ;
        mapObj.writeOpt() = IOobject::AUTO_WRITE;

        if (mapObj.headerOk())
        {
            return new labelIOList
            (
                mapObj
            );
        }
        else
        {
            mapObj.readOpt() = IOobject::NO_READ;

            label defaultRegionI = regionIndex(polyMesh::defaultRegion);

            if (!cellProcMap_[defaultRegionI])
            {
                makeCellProcMap(defaultRegionI);
            }

            if (!cellRegionMap_[regionI])
            {
                makeCellRegionMap(regionI);
            }

            if (!cellProcMap_[regionI])
            {
                makeCellProcMap(regionI);
            }

            const labelIOList& regionMap = *cellRegionMap_[regionI];
            const labelIOList& procMap = *cellProcMap_[regionI];

            labelIOList* cellRegionProcMapPtr =
                new labelIOList
                (
                    mapObj,
                    labelList(procMap.size(),-1)
                );

            labelIOList& cellRegionProcMap = *cellRegionProcMapPtr;

            forAll(procMap,cellI)
            {
                cellRegionProcMap[cellI] =
                    findIndex<labelList>
                    (
                        *cellProcMap_[defaultRegionI],
                        regionMap[procMap[cellI]]
                    );
            }

            return cellRegionProcMapPtr;

//             WarningIn("regionPolyMesh::makeCellRegionProcMap(...)")
//                 << "For the current implementation of the "
//                 << "region-processor cell map, it is NECESSARY, "
//                 << "that the cell-order of cellZones of the mesh of "
//                 << "the default region (region0) of each processor "
//                 << "is equal to the cell-order of the corresponding "
//                 << "sub-region mesh for this processor!"
//                 << endl;
//
//             const label cellZoneRegionIinDefaultI =
//                 mesh(polyMesh::defaultRegion).cellZones().findZoneID
//                 (
//                     regionName(regionI)
//                 );
//
//             if (cellZoneRegionIinDefaultI == -1)
//             {
//                 FatalErrorIn("regionPolyMesh::newCellRegionProcMap(...)")
//                     << "There is no cellZone called "
//                     << regionName(regionI)
//                     << " in mesh of "
//                     << polyMesh::defaultRegion
//                     << abort(FatalError);
//             }
//
//             const labelList& cellZoneRegionIinDefault =
//                 mesh(polyMesh::defaultRegion)
//                 .cellZones()[cellZoneRegionIinDefaultI];
//
//             mapObj.readOpt() = IOobject::READ_IF_PRESENT;
//             mapObj.writeOpt() = IOobject::AUTO_WRITE;
//
//             return new labelIOList
//             (
//                 mapObj,
//                 cellZoneRegionIinDefault
//             );
        }
    }
}

labelIOList* regionPolyMesh::newPointRegionProcMap
(
    const label& regionI
) const
{
    IOobject mapObj
    (
        "pointRegionProcAddressing",
        mesh(regionI).time().constant(),
        polyMesh::meshSubDir,
        mesh(regionI),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    );

    if (regionName(regionI) == polyMesh::defaultRegion)
    {
        return new labelIOList
        (
            mapObj,
            labelList()
        );
    }
    else
    {
        mapObj.readOpt() = IOobject::MUST_READ;
        mapObj.writeOpt() = IOobject::AUTO_WRITE;

        if (mapObj.headerOk())
        {
            return new labelIOList
            (
                mapObj
            );
        }
        else
        {
            mapObj.readOpt() = IOobject::NO_READ;

            label defaultRegionI = regionIndex(polyMesh::defaultRegion);

            if (!pointProcMap_[defaultRegionI])
            {
                makePointProcMap(defaultRegionI);
            }

            if (!pointRegionMap_[regionI])
            {
                makePointRegionMap(regionI);
            }

            if (!pointProcMap_[regionI])
            {
                makePointProcMap(regionI);
            }

            const labelIOList& regionMap = *pointRegionMap_[regionI];
            const labelIOList& procMap = *pointProcMap_[regionI];

            labelIOList* pointRegionProcMapPtr =
                new labelIOList
                (
                    mapObj,
                    labelList(procMap.size(),-1)
                );

            labelIOList& pointRegionProcMap = *pointRegionProcMapPtr;

            forAll(procMap,cellI)
            {
                pointRegionProcMap[cellI] =
                    findIndex<labelList>
                    (
                        *pointProcMap_[defaultRegionI],
                        regionMap[procMap[cellI]]
                    );
            }

            return pointRegionProcMapPtr;
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const labelIOList&
regionPolyMesh::cellMap(const label& regionI) const
{
    if (!cellMap_[regionI])
    {
        makeCellMap(regionI);
    }

    return *cellMap_[regionI];
}

const labelIOList&
regionPolyMesh::cellMap(const word& regionName) const
{
    label regionI = regionIndex(regionName);

    if (!cellMap_[regionI])
    {
        makeCellMap(regionI);
    }

    return *cellMap_[regionI];
}

const labelIOList&
regionPolyMesh::pointMap(const label& regionI) const
{
    if (!pointMap_[regionI])
    {
        makePointMap(regionI);
    }

    return *pointMap_[regionI];
}

const labelIOList&
regionPolyMesh::pointMap(const word& regionName) const
{
    label regionI = regionIndex(regionName);

    if (!pointMap_[regionI])
    {
        makePointMap(regionI);
    }

    return *pointMap_[regionI];
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

