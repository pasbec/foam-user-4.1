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
    solverManagerCoreBase

Description
    ...

\*---------------------------------------------------------------------------*/

#ifndef solverManagerCoreBase_H
#define solverManagerCoreBase_H

#include "solverManagerCoreMesh.H"

// TODO: Implement Storage Items

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template <class MESH> class solverManagerCore;

/*---------------------------------------------------------------------------*\
                         Class solverManagerCore Declaration
\*---------------------------------------------------------------------------*/

template <class MESH>
class solverManagerCoreBase
{
    // Friends
    friend class solverManagerCore<MESH>;

public:

    // Public typedefs

        typedef solverManagerCore<MESH> Manager;

        typedef MESH ManagerMesh;
        typedef typename solverManagerCoreMesh<MESH>::Mesh Mesh;

    // Forward declarations

        class Settings;
        class Storage;
        class Region;
        class Regions;


private:

    // Private data types

        //- Read
        class Read
        {
            //- Friends
            friend class solverManagerCore<MESH>;
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


private:

    // Private data types


        //- Init
        class Init
        {
            //- Friends
            friend class solverManagerCore<MESH>;
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

    // Public data types

        //- Settings
        class Settings : public Read, public Init
        {
            //- Friends
            friend class solverManagerCore<MESH>;
            friend class Region;

        private:

            // Private Member functions

                //- Init sequence
                virtual void init() const {};


        public:

            // Constructor

                Settings() : Read(), Init() {};

            // Destructor

                virtual ~Settings() {};
        };


public:

    // Public data types


        //- Item
        class StorageItem : public Init
        {
            //- Friends
            friend class solverManagerCore<MESH>;
            friend class Storage;

        private:

            // Private Member functions

                //- Read sequence
                virtual void read() const {};


        public:

            // Constructor

                StorageItem() : Read(), Init() {};

            // Destructor

                virtual ~StorageItem() {};
        };


public:

    // Public data types


        //- Storage
        class Storage : public Read, public Init
        {
            //- Friends
            friend class solverManagerCore<MESH>;
            friend class Region;

        private:

            // Private Member functions

                //- Read sequence
                virtual void read() const {};


        public:

            // Constructor

                Storage() : Read(), Init() {};

            // Destructor

                virtual ~Storage() {};
        };


public:

    // Public data types


        //- Region
        class Region : public Read, public Init
        {
            //- Friends
            friend class Regions;

        private:

            // Private Member functions

                //- Read sequence
                virtual void read() const
                {
                    settings().read();
                };

                //- Init sequence
                virtual void init() const
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

                //- Settings
                virtual Settings& settings() const = 0;

                //- Storage
                virtual Storage& storage() const = 0;

                // Basic

                    //- Mesh
                    virtual Mesh& mesh() const = 0;
        };


public:

    // Public data types


        //- Regions
        class Regions : public Read, public Init
        {
            //- Friends
            friend class solverManagerCore<MESH>;

        public:

            // Public typedefs

                typedef HashTable<Region*> RegionHashTable;


        private:

            // Private Data

                //- Hash table
                mutable autoPtr<RegionHashTable> regionHashTable_;

            // Private Member Functions

                //- Read all regions in hash table
                virtual void read() const
                {
                    forAllConstIter(typename RegionHashTable, regionHashTable(), iter)
                    {
                        iter()->read();
                    }
                };

                //- Init all regions in hash table
                virtual void init() const
                {
                    read();

                    forAllConstIter(typename RegionHashTable, regionHashTable(), iter)
                    {
                        iter()->checkInit();
                    }
                };


        protected:

                //- Hash table
                inline RegionHashTable& regionHashTable() const
                {
                    if (regionHashTable_.empty())
                    {
                        regionHashTable_.set(new RegionHashTable());
                    }

                    return regionHashTable_();
                };


        public:

            // Constructor

                Regions() : Read(), Init() {};

            // Destructor

                virtual ~Regions() {};

            // Member Functions

                //- Default region
                virtual Region& defaultRegion() const = 0;

            // Member Operators

                //- Find and return a hashedEntry
                inline Region& operator[](const word& name)
                {
                    return *regionHashTable()[name];
                };
        };
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
