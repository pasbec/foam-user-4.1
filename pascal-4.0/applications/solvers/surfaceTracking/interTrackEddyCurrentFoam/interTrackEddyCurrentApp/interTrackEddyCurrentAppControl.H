/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

Description
    ...

SourceFiles
    interTrackEddyCurrentAppControlI.H
    interTrackEddyCurrentAppControl.C

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#ifndef interTrackEddyCurrentAppControl_H
#define interTrackEddyCurrentAppControl_H

#include "interTrackEddyCurrentAppInclude.H"

#include "solutionControl.H"
#include "solverPerformanceData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace interTrackEddyCurrentApp
{

/*---------------------------------------------------------------------------*\
                   Class interTrackEddyCurrentApp::Control Declaration
\*---------------------------------------------------------------------------*/

class Control
:
    public solutionControl
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        Control(const Control&);

        //- Disallow default bitwise assignment
        void operator=(const Control&);


protected:

    // Protected structs

        //- Electromagnetic settings
        struct emSettings
        {
            Switch enabled;
            int outputTimeIndexCycle;
            int timeIndexCycle;
            scalar timeCycle;
            scalar relDeltaAmax;
        };

        //- Electromagnetic update data
        struct emUpdateData
        {
            Switch update;
            int counter;
            int outputTimeIndex;
            scalar lastTime;
        };

    // Protected data

        //- Mesh
        const regionDynamicFvMesh& mesh_;

        // Electromagnetic

            //- Previous cell centres
            volVectorField& prevC_;

            //- Estimated relative deviation of A
            volScalarField& relDeltaA_;

            //- Properties dictionary
            const dictionary& propertiesDict_;

            //- Electromagnetic settings dictionary
            const dictionary& emUpdateSettingsDict_;

            //- Electromagnetic settings
            emSettings emSettings_;

            //- Electromagnetic update dictionary
            IOdictionary emUpdateDataDict_;

            //- Electromagnetic update data
            emUpdateData emUpdateData_;

    // Protected Member Functions

        //- Read controls
        virtual void read();

        //- Return true if an electromagnetic update is not necessary
        virtual bool criteriaSatisfied();

        // Update checks

            //- Counter is zero
            inline bool updateZeroCounter() const;

            //- Output time cycle interval is due
            inline bool updateOutputTimeIndex() const;

            //- Time index cycle interval is due
            inline bool updateTimeIndex() const;

            //- Time cycle interval is due
            inline bool updateTime() const;

            //- Relative A deviation has grown to large
            bool updateRelDeltaA() const;

            //- Update is due
            bool update();


public:

    //- Run-time type information
    TypeName("interTrackEddyCurrentApp::Control");


    // Constructors

        //- Construct from mesh and the name of control sub-dictionary
        Control
        (
            regionDynamicFvMesh& mesh,
            const dictionary& propertiesDict,
            const word& dictName="INTERTRACKEDDYCURRENT"
        );


    //- Destructor
    virtual ~Control();


    // Member Functions

        // Evolution

            //- Main control loop
            virtual bool loop();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interTrackEddyCurrentApp


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interTrackEddyCurrentAppControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "interTrackEddyCurrentAppControl.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
