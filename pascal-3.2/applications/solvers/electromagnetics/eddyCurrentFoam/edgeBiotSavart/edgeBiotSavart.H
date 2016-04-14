/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

Class
    Foam::edgeBiotSavart

Description
    ...

SourceFiles
    edgeBiotSavart.C

\*---------------------------------------------------------------------------*/

#ifndef edgeBiotSavart_H
#define edgeBiotSavart_H

#include "IOdictionary.H"
#include "fvMesh.H"
#include "volFields.H"
#include "featureEdgeMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class edgeBiotSavart Declaration
\*---------------------------------------------------------------------------*/

class edgeBiotSavart
:
    public IOdictionary
{
public:

    //- Real/Imaginary part enums
    enum ComplexPart
    {
        REAL = 0,
        IMAG = 1
    };

    //- Inductor data
    struct InductorData
    {
        word name;
        fileName file;
        Switch reverse;
        scalar current;
        scalar phase;
    };


private:

    // Private data

        //- Mesh
        const fvMesh& mesh_;

        //- List of inductors
        List<InductorData> inductorData_;

        //- List of inductor meshes
        List<featureEdgeMesh*> inductorMeshes_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        edgeBiotSavart(const edgeBiotSavart&);

        //- Disallow default bitwise assignment
        void operator=(const edgeBiotSavart&);

        // Core Calculation

            // Magnetic vector potential

                //- Calculate A for given points
                //  Parameter part may be used to respect
                //  corresponding phase angle
                tmp<vectorField> calcA
                (
                    const pointField& points,
                    ComplexPart part = REAL
                ) const;

            // Magnetic field density

                //- Calculate B for given points
                tmp<vectorField> calcB
                (
                    const pointField& points,
                    ComplexPart part = REAL
                ) const;


public:

    // Runtime type information
    TypeName("edgeBiotSavart");

    // Constructor

        edgeBiotSavart(const fvMesh& mesh);

    // Destructor

        virtual ~edgeBiotSavart();

    // Member Functions

        // Magnetic vector potential

            //- Calculate A for internal cell centres
            //  Parameter part may be used to respect
            //  corresponding phase angle
            tmp<vectorField> internalA
            (
                ComplexPart part = REAL
            ) const;

            //- Calculate A for boundary patch faces centres
            //  Parameter part may be used to respect
            //  corresponding phase angle
            tmp<vectorField> boundaryPatchA
            (
                const label& patchI,
                ComplexPart part = REAL
            ) const;

            //- Calculate A for boundary patch faces centres
            //  Parameter part may be used to respect
            //  corresponding phase angle
            tmp<vectorField> boundaryPatchA
            (
                const word& patchName,
                ComplexPart part = REAL
            ) const;

            //- Calculate A
            //  Parameter part may be used to respect
            //  corresponding phase angle
            tmp<volVectorField> A
            (
                ComplexPart part = REAL
            ) const;

        // Magnetic field density

            //- Calculate B for internal cell centres
            //  Parameter part may be used to respect
            //  corresponding phase angle
            tmp<vectorField> internalB
            (
                ComplexPart part = REAL
            ) const;

            //- Calculate B for boundary patch faces centres
            //  Parameter part may be used to respect
            //  corresponding phase angle
            tmp<vectorField> boundaryPatchB
            (
                const label& patchI,
                ComplexPart part = REAL
            ) const;

            //- Calculate B for boundary patch faces centres
            //  Parameter part may be used to respect
            //  corresponding phase angle
            tmp<vectorField> boundaryPatchB
            (
                const word& patchName,
                ComplexPart part = REAL
            ) const;

            //- Calculate B
            //  Parameter part may be used to respect
            //  corresponding phase angle
            tmp<volVectorField> B
            (
                ComplexPart part = REAL
            ) const;

        //- Write inductor/mesh data to vtk
        void writeVTK() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
