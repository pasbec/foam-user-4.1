/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    meshFluxesObject

Description
    FunctionObject for mesh fluxes.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    meshFluxesObject.C

\*---------------------------------------------------------------------------*/

#ifndef meshFluxesObject_H
#define meshFluxesObject_H

#include "functionObject.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class meshFluxesObject Declaration
\*---------------------------------------------------------------------------*/

class meshFluxesObject
:
    public functionObject
{
    // Private data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time& time_;

        //- Region name
        word regionName_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        meshFluxesObject(const meshFluxesObject&);

        //- Disallow default bitwise assignment
        void operator=(const meshFluxesObject&);

        //- Triangle swept volume
        scalar triSweptVol
        (
            const triPointRef& oT,
            const triPointRef& nT,
            const scalar deltaT,
            const label method
        ) const;

        //- Compute face swept-volume by central decomposition
        scalar sweptVol
        (
            const face& f,
            const pointField& oldPoints,
            const pointField& newPoints,
            const scalar deltaT,
            const label method
        ) const;

        //- Compute mesh-fluxes by method
        tmp<scalarField> sweptVols
        (
            const fvMesh& mesh,
            const pointField& newPoints,
            const pointField& oldPoints,
            const scalar deltaT,
            const label method
        ) const;

        // Set fluxes into GeometricField
        void setFluxes
        (
            const scalar& rDeltaT,
            const scalarField& sweptVols,
            surfaceScalarField& phi
        ) const;

public:

    //- Runtime type information
    TypeName("meshFluxes");


    // Constructors

        //- Construct from components
        meshFluxesObject
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Member Functions

        //- start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop
        virtual bool execute();

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
