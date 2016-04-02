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
    solverManager

Description
    ...

SourceFiles
    solverManager.C

\*---------------------------------------------------------------------------*/

#ifndef solverManager_H
#define solverManager_H

#include "autoPtr.H"
#include "argList.H"
#include "foamTime.H"
#include "IOdictionary.H"
#include "solutionControl.H"

#include "solverManagerMessages.H"
#include "solverManagerBaseNested.H"
#include "solverManagerRegions.H"
#include "solverManagerRegion.H"
#include "solverManagerStorage.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// TODO: Remake setManagerScope, ...


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class solverManager Declaration
\*---------------------------------------------------------------------------*/

template <class MESH>
class solverManager
{
public:

    // Public typedefs

        typedef solverManager Manager;
        typedef MESH ManagerMesh;
        typedef typename MESH::Mesh ThisMesh;

        typedef solverManagerMessages<MESH> Messages;
        typedef solutionControl Control;

    // Public Nested Classes

        class Regions {};
        class Region {};
        class Storage {};
        class Item {};

    // Public Structs

// TODO: Add global settings

// TODO: Rename parameters
        //- Time parameters
        struct parameters
        {
            bool adjustTimeStep;
            scalar maxCo;
            scalar maxDeltaT;
            scalar CoNum;
        };


protected:

    // Protected Member Functions

        // Evolution

            //- Set time step width
            //  A return value of true means active
            virtual bool setDeltaT(scalar& deltaT) const
            {
                return false;
            };

            //- Set time parameters
            //  A return value of true means active
            virtual bool setParameters(parameters& tp) const
            {
                return false;
            };

            //- Set Courant Number
            //  A return value of true means active
            virtual bool setCoNum(scalar& CoNum) const
            {
                return false;
            };

            //- Init sequence
            virtual void init() const = 0;

            //- Read sequence
            virtual void read() const = 0;

            //- Finalize sequence
            virtual void finalize() const = 0;

            //- Pre write sequence
            virtual void writePre() const {};

            //- Post write sequence
            virtual void writePost() const {};


private:

    // Private Data

        //- Time parameters
        mutable autoPtr<parameters> paramPtr_;

        //- Messages
        mutable autoPtr<Messages> msgPtr_;

        // Basic

            //- Reference to args
            const argList& args_;

            //- Reference to time
            Time& time_;

            //- Reference to mesh
            ManagerMesh& mesh_;

            //- Master manager?
            const bool master_;

            //- Properties dictionary
            IOdictionary propertiesDict_;

        // Evolution

            //- Pre runTime phase?
            mutable bool prePhase_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        solverManager(const solverManager&);

        //- Disallow default bitwise assignment
        void operator=(const solverManager&);

        // Basic

            //- Error if not master manager
            void errorIfNotMaster() const;

        // Evolution

            //- Pre runTime phase?
            inline bool& prePhase() const
            {
                return prePhase_;
            };

            //- Read time parameters
            void readParameters(parameters& tc) const;

            //- Time step setting
            void calcDeltaT(scalar& deltaT) const;

            //- Time step application
            void applyDeltaT() const;


public:

    // Constructor

        solverManager
        (
            const argList& args,
            Time& time,
            ManagerMesh& mesh,
            const word& prefix,
            const bool& master = true
        );

    // Destructor

        virtual ~solverManager() {};

    // Member Functions

        //- Time parameters
        const parameters& param() const;

        //- Messages
        const Messages& msg() const;

        //- Regions
        virtual Regions& regions() const = 0;

        // Basic

            //- Reference to args
            inline const argList& args() const
            {
                return args_;
            };

            //- Reference to time
            inline Time& time() const
            {
                return time_;
            };

            //- Reference to mesh
            inline ManagerMesh& mesh() const
            {
                return mesh_;
            };

            //- Is this the master manager?
            inline const bool& master() const
            {
                return master_;
            };

            //- Refrence to control dictionary
            inline const dictionary& controlDict() const
            {
                return time_.controlDict();
            }

            //- Reference to properties dictionary
            inline const dictionary& propertiesDict() const
            {
                return propertiesDict_;
            };

            //- Reference to settings dictionary
            inline const dictionary& settingsDict() const
            {
                return propertiesDict_.subDict("settings");
            };

            //- Reference to settings dictionary for named region
            const dictionary& settingsRegionDict
            (
                const word& region = polyMesh::defaultRegion
            ) const;

            //- Reference to storage dictionary
            inline const dictionary& storageDict() const
            {
                return propertiesDict_.subDict("storage");
            };

            //- Reference to storage dictionary for named region
            const dictionary& storageRegionDict
            (
                const word& region = polyMesh::defaultRegion
            ) const;

            //- Reference to dictionary of named item of storage
            bool storageItemDict
            (
                const word& item,
                dictionary& dict,
                const word& region = polyMesh::defaultRegion
            ) const;

        // Evolution

            //- Main time loop wrapper
            bool loop() const;

            //- Main time run wrapper
            bool run() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "solverManager.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
