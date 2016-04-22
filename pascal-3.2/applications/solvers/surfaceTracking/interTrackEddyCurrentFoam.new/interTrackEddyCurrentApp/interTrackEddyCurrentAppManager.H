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
    Foam::interTrackEddyCurrentApp::Manager

Description
    ...

SourceFiles
    interTrackEddyCurrentAppManager.C

\*---------------------------------------------------------------------------*/

#ifndef interTrackEddyCurrentAppManager_H
#define interTrackEddyCurrentAppManager_H

#include "interTrackEddyCurrentAppInclude.H"

#include "interTrackAppManager.H"
#include "eddyCurrentAppManager.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace interTrackEddyCurrentApp
{

/*---------------------------------------------------------------------------*\
               Class interTrackEddyCurrentApp::Manager Declaration
\*---------------------------------------------------------------------------*/

class Manager
:
    public solverManager<regionDynamicFvMesh>
{
public:

    // Public typedefs and forward declarations

        typedef interTrackApp::Control Control;
        class Settings;
        class Storage;

    // Solver Manager Construction

        SM_GLOBALSETTINGS_NEW();
        SM_GLOBALSETTINGS_ADD();


        SM_GLOBALSTORAGE_NEW();
        SM_GLOBALSTORAGE_ADD();


        SM_REGIONS_NEW();
        SM_REGIONS_ADD();


private:

    // Private data

        //- Pointer for manager of interTrackApp
        mutable autoPtr<interTrackApp::Manager> interTrackAppManager_;

        //- Pointer for manager of eddyCurrentApp
        mutable autoPtr<eddyCurrentApp::Manager> eddyCurrentAppManager_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        Manager(const Manager&);

        //- Disallow default bitwise assignment
        void operator=(const Manager&);


public:

    // Runtime type information
    TypeName("interTrackEddyCurrentApp::Manager");

    // Constructor

        Manager
        (
            const argList& args,
            Time& time,
            regionDynamicFvMesh& mesh,
            bool master = true,
            const word& prefix = "interTrackEddyCurrent"
        );

    // Destructor

        virtual ~Manager();

    // Public Member Functions

        //- Control
        virtual Control& control() const;

        // Managers

            //- interTrackApp
            interTrackApp::Manager& interTrackAppManager() const;

            //- eddyCurrentApp
            eddyCurrentApp::Manager& eddyCurrentAppManager() const;

        // Evolution

            //- Set Courant Number
            virtual bool setCoNum(scalar& CoNum) const;

            //- Init sequence
            virtual void init() const;

            //- Read sequence
            virtual void read() const;

            //- Next sequence
            virtual void next() const;

            //- Write sequence
            virtual void write() const;

            //- Finalize sequence
            virtual void finalize() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interTrackEddyCurrentApp

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "interTrackEddyCurrentAppManager.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
