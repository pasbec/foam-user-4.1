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
    Foam::solutionData

Description
    Database for solution data, solver performance and other reduced data.

SourceFiles
    solutionData.C
    solutionDataTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef solutionData_H
#define solutionData_H

#include "IOdictionary.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class solutionData Declaration
\*---------------------------------------------------------------------------*/

class solutionData
:
    public IOdictionary
{
    // Private data

        //- Previously used time-index, used for reset between iterations
        mutable label prevTimeIndex_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        solutionData(const solutionData&);

        //- Disallow default bitwise assignment
        void operator=(const solutionData&);


public:

    //- Runtime type information
    TypeName("solutionData");


    // Constructors

        //- Construct from mesh
        solutionData(const fvMesh& mesh);


    // Member Functions

        // Access

            //- Return the dictionary of solver performance data
            //  which includes initial and final residuals for convergence
            //  checking
            const dictionary& solverPerformanceDict() const;

            //- Add/set the solverPerformance entry for the named field
            template<class spType>
            void setSolverPerformance
            (
                const word& name,
                const spType& sp
            ) const;

            //- Add/set the solverPerformance entry using its name
            template<class spType>
            void setSolverPerformance
            (
                const spType& sp
            ) const;

// TEST TODO: Move/Wrap these templates to/in solutionControl
            //- Solve tmp<matrix> and store solverPerformance in database
            template<class matrixType>
            void solve
            (
                const tmp<matrixType>& tm
            ) const;

// TEST TODO: Move/Wrap these templates to/in solutionControl
            //- Solve matrix and store solverPerformance in database
            template<class matrixType>
            void solve
            (
                matrixType& m
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "solutionDataTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
