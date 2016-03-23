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

Class
    eddyCurrentControl

Description
    ...

SourceFiles
    eddyCurrentControl.C

\*---------------------------------------------------------------------------*/

#ifndef eddyCurrentControl_H
#define eddyCurrentControl_H

#include "PstreamReduceOps.H"
#include "mathematicalConstants.H"
#include "physicalConstants.H"
#include "solverControl.H"
#include "regionModelling.H"

#include "fixedGradientFvPatchFields.H"
#include "tangentialMagneticFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class eddyCurrentControl Declaration
\*---------------------------------------------------------------------------*/

class eddyCurrentControl
:
    public solverControl<regionFvMesh>
{
private:

    // Private Data

        //- Storage

        //- Controls

        //- Properties

            //- Name of conducting region
            word conductorRegionName_;

            //- Label of conducting region
            label conductorRegion_;

            //- Name of conductor interface
            word interfacePatchName_;

            //- Label of conductor interface
            word interfacePatch_;

            //- Frequency
            dimensionedScalar frequency_;

            //- Cyclic velocity based on frequency
            dimensionedScalar omega_;

        //- Dictionaries

            //-Solution dictionary of conducting region
            solution& conductorSolutionDict_;

            //- AV-solution dictionary
            dictionary AVdict_;

            //- A-solution dictionary
            dictionary Adict_;

            //- V-solution dictionary
            dictionary Vdict_;

        //- Solution properties

            //- Solution tolerance
            scalar tol_;

            //- Relative solution tolerance
            scalar relTol_;

            //- Outer maximum number of iterations
            int maxIter_;

            //- Realxation factor
            scalar relax_;

            //- Is mesh 3D?
            mutable bool mesh3D_;

            //- Solution direction
            mutable Vector<label> solDir_;

        //- Solution state variables

            //- Stop state
            mutable bool stop_;

            //- Loop/Convergence state
            mutable bool loop_;

            //- Iteration counter
            mutable int iter_;

            //- Sub-loop/Convergence state
            mutable bool subLoop_;

            //- Sub-iteration counter
            mutable int subIter_;

            //- Sub-iteration relative tolerance
            mutable scalar subRelTol_;

            //- Sub-iteration outer relative solution
            //-  tolerance limit for subRelTol down-scaling
            mutable scalar subRelTolLim_;

            //- Sub-iteration relative tolerance
            mutable scalar subScale_;

            //- Sub-iteration relax
            mutable scalar subRelax_;

            //- Temporary copy of A-solution dictionary
            mutable dictionary subAdict_;

            //- Temporary copy of V-solution dictionary
            mutable dictionary subVdict_;

        //- Solution residual variables

            //- AV-residual
            mutable scalar AVres_;

            //- Relative residual
            mutable scalar AVrelRes_;

            //- A-residual
            mutable scalar Ares_;

            //- V-residual
            mutable scalar Vres_;

            //- Old AV-residual
            mutable scalar oldAVres_;

            //- Old relative residual
            mutable scalar oldAVrelRes_;

            //- Old AV-residual switch
            mutable bool oldAVresSet_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        eddyCurrentControl(const eddyCurrentControl&);

        //- Disallow default bitwise assignment
        void operator=(const eddyCurrentControl&);

        //- Inline helper

            //- Check if relax factor is valid
            inline bool checkRelax() const
            {
                return (relax_ > 0.0) && (relax_ < 1.0);
            };

            //- Check if relax factor is valid
            inline bool checkSubRelax() const
            {
                return (subRelax_ > 0.0) && (subRelax_ < 1.0);
            };

            //- Iteration only initialized?
            inline bool noIter() const
            {
                return (iter_ == -1);
            };

            //- Is this the first iteration?
            inline bool firstIteration() const
            {
                return (iter_ == 0 && subIter_ == 0);
            };

            //- Is this the first sub-iteration?
            inline bool firstSubIteration() const
            {
                return (subIter_ == 0);
            };

            //- Check minimum iterations state
            inline bool checkMaxIterations() const
            {
                return (iter_ < maxIter_);
            };

            //- Check convergence state
            inline bool checkConvergence() const
            {
                return (AVres_ <= tol_);
            };

            //- Check relative convergence state
            inline bool checkSubConvergence() const
            {
                return (max(oldAVrelRes_, AVrelRes_) <= subRelTolLim_);
            };

        //- Helper

            //- Read dictionaries and update data
            void readDictDataIfModified();

            //- Check convergence state
            void updateResidual() const;

            //- Update sub-disctionaries data
            void resetSubDictionaries() const;

            //- Decrease sub-tolerance
            void decreaseSubTolerance() const;

            //- Decrease sub-scale
            void decreaseSubScale() const;

            //- Check and return sub-loop state
            const bool& subLoop() const;

            //- Reset
            void reset();

public:

    // Runtime type information
    TypeName("eddyCurrentControl");

    // Constructors

        //- Construct from mesh and name
        eddyCurrentControl
        (
            const argList& args,
            Time& time,
            regionFvMesh& mesh,
            const word& name = "eddyCurrent",
            const bool& master = true
        );

    // Destructor

        virtual ~eddyCurrentControl() {};

    // Member Functions

        //- Storage
        virtual storage& data()
        {
            notImplemented(type() + "::data()");
            storage* storagePtr_ = NULL;
            return *storagePtr_;
        };

        //- Controls

        //- Properties

            //- Reference to conductor region name
            const word& conductorName() const
            {
                return conductorRegionName_;
            };

            //- Reference to conductor region ID
            const label& conductor() const
            {
                return conductorRegion_;
            };

            //- Reference to conductor interface name
            const word& interfaceName() const
            {
                return interfacePatchName_;
            };

            //- Reference to conductor interface ID
            const word& interface() const
            {
                return interfacePatch_;
            };

            //- Reference to region control
            const dimensionedScalar& frequency() const
            {
                return frequency_;
            };

            //- Reference to region control
            const dimensionedScalar& omega() const
            {
                return omega_;
            };

        //- Solution properties

            //- Mesh geometry is 3D
            bool meshIs3D() const
            {
                return mesh3D_;
            };

            //- Mesh solution direction
            const Vector<label>& solutionDir() const
            {
                return solDir_;
            };

            //- Reference to solution tolerance
            const scalar& tolerance() const
            {
                return tol_;
            };

        //- Solution state variables

            //- Reference to temporary A-solution dictionary
            const dictionary& Adict() const
            {
                return subAdict_;
            };

            //- Reference to temporary V-solution dictionary
            const dictionary& Vdict() const
            {
                return subVdict_;
            };

            //- Relax V-volVectorField
            void relax(volVectorField& vf) const;

            //- Relax V-regionVolVectorField
            void relax(const regionVolVectorField& rvf) const;

            //- Is this the first iteration?
            bool firstIter() const
            {
                return firstIteration();
            };

            //- Is this the first sub-iteration?
            bool firstSubIter() const
            {
                return firstSubIteration();
            };

            //- Check, set and return loop state
            const bool& loop();

        //- Solution residual variables

            //- Set A-residual
            void setResidualOfA(scalar residual) const;

            //- Set V-residual
            void setResidualOfV(scalar residual) const;

        //- Time and database

            //- WriteData member function required by regIOobject
            bool writeData(Ostream&) const
            {
                return false;
            };
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "eddyCurrentControl.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
