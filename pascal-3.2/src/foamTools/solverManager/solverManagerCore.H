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
    solverManagerCore

Description
    ...

SourceFiles
    solverManagerCore.C

\*---------------------------------------------------------------------------*/

#ifndef solverManagerCore_H
#define solverManagerCore_H

#include "autoPtr.H"
#include "HashPtrTable.H"
#include "argList.H"
#include "foamTime.H"
#include "IOdictionary.H"
#include "solutionControl.H"

#include "solverManagerCoreName.H"
#include "solverManagerCoreMesh.H"
#include "solverManagerCoreBase.H"
#include "solverManagerCoreMessages.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class solverManagerCore Declaration
\*---------------------------------------------------------------------------*/

template <class MESH>
class solverManagerCore
:
    private solverManagerCoreName,
    private solverManagerCoreMessages<MESH>,
    private solverManagerCoreBase<MESH>
{
public:

    // Public typedefs

        typedef solverManagerCore<MESH> Manager;

        typedef MESH ManagerMesh;
        typedef typename solverManagerCoreBase<MESH>::Mesh Mesh;

        typedef solutionControl Control;
        typedef typename solverManagerCoreMessages<MESH>::Messages Messages;
        typedef typename solverManagerCoreBase<MESH>::Settings Settings;
        typedef typename solverManagerCoreBase<MESH>::Storage Storage;
        typedef typename solverManagerCoreBase<MESH>::Region Region;
        typedef typename solverManagerCoreBase<MESH>::Regions Regions;


public:

    // Public data types

        //- Parameters
        class Parameters
        {
        public:

            // Public structs

                //- Time parameters
                struct Time
                {
                    bool adjustTimeStep;
                    scalar maxCo;
                    scalar maxDeltaT;
                    scalar CoNum;
                };

            // Data

                //- Time parameters
                Time time;
        };

    // Public typedefs

        typedef typename Parameters::Time TimeParameters;


private:

    // Private Data

        //- Parameters
        mutable autoPtr<Parameters> parametersPtr_;

        //- Messages
        mutable autoPtr<Messages> messagesPtr_;

        // Basic

            //- Reference to args
            const argList& args_;

            //- Reference to time
            Time& time_;

            //- Reference to mesh
            ManagerMesh& mesh_;

            //- Master manager?
            const bool master_;

        // Dictionaries

            //- Properties dictionary
            IOdictionary propertiesDict_;

        // Evolution

            //- Pre runTime phase?
            mutable bool prePhase_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        solverManagerCore(const solverManagerCore&);

        //- Disallow default bitwise assignment
        void operator=(const solverManagerCore&);

        // Basic

            //- Error if not master manager
            void errorIfNotMaster() const;

        // Evolution

            //- Pre runTime phase?
            inline bool& prePhase() const
            {
                return prePhase_;
            };

            //- Read parameters
            void readParameters(Parameters& par) const;

            //- Read time parameters
            void readTimeParameters(TimeParameters& timePar) const;

            //- Time step setting
            void calcDeltaT(scalar& deltaT) const;

            //- Time step application
            void applyDeltaT() const;


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
            virtual bool setTimeParameters(TimeParameters& timePar) const
            {
                return false;
            };

            //- Set Courant Number
            //  A return value of true means active
            virtual bool setCoNum(scalar& CoNum) const
            {
                return false;
            };

            //- Next sequence
            virtual void next() const = 0;

            //- Write sequence
            virtual void write() const = 0;

            //- Finalize sequence
            virtual void finalize() const = 0;


public:

    // Constructor

        solverManagerCore
        (
            const argList& args,
            Time& time,
            ManagerMesh& mesh,
            const word& prefix,
            const bool& master = true
        );

    // Destructor

        virtual ~solverManagerCore() {};

    // Member Functions

        //- Parameters
        const Parameters& parameters() const;

        //- Messages
        const Messages& messages() const;

        //- Control
        virtual Control& control() const = 0;

        //- Settings
        virtual Settings& settings() const = 0;

        //- Storage
        virtual Storage& storage() const = 0;

        //- Regions
        virtual Regions& regions() const = 0;

        // Basic

            //- Main arguments
            inline const argList& args() const
            {
                return args_;
            };

            //- Time
            inline Time& time() const
            {
                return time_;
            };

            //- Mesh
            inline ManagerMesh& mesh() const
            {
                return mesh_;
            };

            //- Is this the master manager?
            inline const bool& master() const
            {
                return master_;
            };

        // Dictionaries

            //- Control dictionary
            inline const dictionary& controlDict() const
            {
                return time_.controlDict();
            }

            //- Properties dictionary
            inline const dictionary& propertiesDict() const
            {
                return propertiesDict_;
            };

            inline const dictionary& settingsDict() const
            {
                return propertiesDict_.subDict("settings");
            };

            //- Storage dictionary
            inline const dictionary& storageDict() const
            {
                return propertiesDict_.subDict("storage");
            };


        // Evolution

            //- Read sequence
            void read() const;

            //- Init sequence
            void init() const;

            //- Main time once wrapper
            bool once() const;

            //- Main time loop wrapper
            bool loop() const;

            //- Main time run wrapper
            bool run() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "solverManagerCore.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
