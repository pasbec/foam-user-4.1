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

\*---------------------------------------------------------------------------*/

#include "solutionData.H"
#include "foamTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class spType>
void Foam::solutionData::setSolverPerformance
(
    const word& name,
    const spType& sp
) const
{
    dictionary& dict = const_cast<dictionary&>(solverPerformanceDict());

    List<dictionary> fieldSpDicts;

    if (prevTimeIndex_ != this->time().timeIndex())
    {
        // reset solver performance between iterations
        prevTimeIndex_ = this->time().timeIndex();
        dict.clear();
    }
    else
    {
        dict.readIfPresent(name, fieldSpDicts);
    }

    dictionary newSpDict;
    {
        newSpDict.add("solverName", sp.solverName());
        newSpDict.add("name", name);
        newSpDict.add("initialResidual", sp.initialResidual());
        newSpDict.add("finalResidual", sp.finalResidual());
        newSpDict.add("nIterations", sp.nIterations());
        newSpDict.add("converged", sp.converged());
        newSpDict.add("singular", sp.singular());
    }

    fieldSpDicts.setSize
    (
        fieldSpDicts.size()+1,
        newSpDict
    );

    dict.set(name, fieldSpDicts);
}


template<class spType>
void Foam::solutionData::setSolverPerformance
(
    const spType& sp
) const
{
    setSolverPerformance(sp.name(), sp);
}


// TEST TODO: Move/Wrap these templates to/in solutionControl
template<class matrixType>
void Foam::solutionData::solve
(
    const tmp<matrixType>& tm
) const
{
    matrixType& m = *tm.ptr();

    setSolverPerformance
    (
        m.psi().name(),
        m.solve()
    );

    tm.clear();
}


// TEST TODO: Move/Wrap these templates to/in solutionControl
template<class matrixType>
void Foam::solutionData::solve
(
    matrixType& m
) const
{
    setSolverPerformance
    (
        m.psi().name(),
        m.solve()
    );
}


// ************************************************************************* //
