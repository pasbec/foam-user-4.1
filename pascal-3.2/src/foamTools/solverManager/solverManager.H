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
#include "HashPtrTable.H"
#include "argList.H"
#include "foamTime.H"
#include "IOdictionary.H"
#include "dynamicFvMesh.H"
#include "solutionControl.H"

#include "solverManagerMessages.H"
#include "solverManagerWrapper.H"
#include "solverManagerGlobalSettings.H"
#include "solverManagerSettings.H"
#include "solverManagerStorage.H"
#include "solverManagerRegion.H"
#include "solverManagerRegions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// TODO: Remake setManagerScope, ...

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class solverManagerMesh Declaration
\*---------------------------------------------------------------------------*/

template <typename MESH>
struct solverManagerMesh
{
    typedef typename MESH::Mesh Mesh;
};

template <>
struct solverManagerMesh<dynamicFvMesh>
{
    typedef dynamicFvMesh Mesh;
};

// TODO: Specialization for regionDynamicFvMesh

/*---------------------------------------------------------------------------*\
                         Class solverManager Declaration
\*---------------------------------------------------------------------------*/

template <class MESH>
class solverManager
{
public:

    // Public typedefs

        typedef solverManager<MESH> Manager;
        typedef MESH ManagerMesh;
        typedef typename solverManagerMesh<MESH>::Mesh Mesh;

        typedef solverManagerMessages<MESH> Messages;
        typedef solutionControl Control;

    //- Forward declarations

        class Settings;
        class Storage;
        class Region;
        class Regions;


private:

    // Private Nested Classes

        //- Read
        class Read
        {
            friend class solverManager<MESH>;
            friend class Settings;
            friend class Storage;
            friend class Region;
            friend class Regions;


        private:

            // Private Member functions

                //- Init sequence
                virtual void read() const = 0;


        public:

            // Constructor

                Read() {};

            // Destructor

                virtual ~Read() {};
        };


        //- Init
        class Init
        {
            friend class Settings;
            friend class Storage;
            friend class Region;
            friend class Regions;


        private:

            // Private Data

                //- Initalized?
                const bool initialized_;

            // Private Member functions

                //- Init sequence
                virtual void init() const = 0;


        protected:

            // Protected Member Functions

                //- Check init wrapper
                bool checkInit() const
                {
                    if (!initialized_)
                    {
                        init();

                        bool& initialized =
                            const_cast<bool&>(initialized_);
                        initialized = true;
                    }

                    return initialized_;
                }


        public:

            // Constructor

                Init() : initialized_(false) {};

            // Destructor

                virtual ~Init() {};

            // Member Functions

                //- Initialized?
                inline bool initialized() const
                {
                    return initialized_;
                }
        };


public:

    // Public Nested Classes

        //- Settings
        class Settings : public Read, public Init
        {
        private:

            // Private Member functions

                //- Init sequence
                void init() const {};


        public:

            // Constructor

                Settings() : Read(), Init() {};

            // Destructor

                virtual ~Settings() {};
        };


// TODO: Items
//         //- Item
//         class Item
//         {
//         };


        //- Storage
        class Storage : public Read, public Init
        {
        public:

            // Public typedefs

// TODO: Items
//                 typedef typename solverManager<MESH>::Item Item;


        private:

            // Private Member functions

                //- Read sequence
                void read() const {};


        public:

            // Constructor

                Storage() : Read(), Init() {};

            // Destructor

                virtual ~Storage() {};
        };


        //- Region
        class Region : public Read, public Init
        {
            //- Friends
            friend class Regions;

        public:

            // Public typedefs

                typedef typename solverManager<MESH>::Control Control;
                typedef typename solverManager<MESH>::Settings Settings;
                typedef typename solverManager<MESH>::Storage Storage;


        private:

            // Private Member functions

                //- Read sequence
                void read() const
                {
                    settings().read();
                };

                //- Init sequence
                void init() const
                {
                    read();
                    storage().checkInit();
                };


        public:

            // Constructor

                Region() : Read(), Init()  {};

            // Destructor

                virtual ~Region() {};

            // Member Functions

                //- Control
                virtual Control& control() const = 0;

                //- Settings
                virtual Settings& settings() const = 0;

                //- Storage
                virtual Storage& storage() const = 0;
        };


        //- Regions
        class Regions : public Read, public Init
        {
            //- Friends
            friend class solverManager<MESH>;


        public:

            // Public typedefs

                typedef typename solverManager<MESH>::Region Region;
                typedef HashPtrTable<Region> RegionHashTable;


        private:

            // Private Data

                //- Hash table
                mutable autoPtr<RegionHashTable> regions_;

            // Private Member Functions

                //- Read all regions in hash table
                virtual void read() const
                {
                    forAllConstIter(typename RegionHashTable, regions(), iter)
                    {
                        iter()->read();
                    }
                };

                //- Init all regions in hash table
                virtual void init() const
                {
                    read();

                    forAllConstIter(typename RegionHashTable, regions(), iter)
                    {
                        iter()->checkInit();
                    }
                };


        protected:

                //- Hash table
                inline RegionHashTable& regions() const
                {
                    if (regions_.empty())
                    {
                        regions_.set(new RegionHashTable());
                    }

                    return regions_();
                };


        public:

            // Constructor

                Regions() : Read(), Init() {};

            // Destructor

                virtual ~Regions() {};

            // Member Functions

                //- Default region
                virtual Region& defaultRegion() const = 0;
        };


public:

    // Private Nested Classes

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

        typedef Settings GlobalSettings;
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

        //- Parameters
        const Parameters& parameters() const;

        //- Messages
        const Messages& messages() const;

        //- Global Settings
        virtual GlobalSettings& globalSettings() const = 0;

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

            inline const dictionary& settingsDict() const
            {
                return propertiesDict_.subDict("settings");
            };

            //- Reference to storage dictionary
            inline const dictionary& storageDict() const
            {
                return propertiesDict_.subDict("storage");
            };


        // Evolution

            //- Read sequence
            void read() const;

            //- Init sequence
            void init() const;

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
