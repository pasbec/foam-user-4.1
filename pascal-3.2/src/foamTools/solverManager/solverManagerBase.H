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
    solverManagerBase

Description
    ...

\*---------------------------------------------------------------------------*/

#ifndef solverManagerBase_H
#define solverManagerBase_H

#include "dynamicFvMesh.H"
#include "regionDynamicFvMesh.H"

#include "solutionControl.H"

// TODO: Implement Storage Items

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template <class MESH> class solverManager;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solverManagerBase
{

template <class MESH> class MeshType;
template <class MESH> class Settings;
template <class MESH> class Storage;
template <class MESH> class Region;
template <class MANAGERMESH, class MESH> class Regions;

/*---------------------------------------------------------------------------*\
                  Class solverManagerBase::Mesh Declaration
\*---------------------------------------------------------------------------*/

template <typename MESH>
struct MeshType
{
    typedef typename MESH::Mesh Type;
};

template <>
struct MeshType<dynamicFvMesh>
{
    typedef dynamicFvMesh Type;
};

template <>
struct MeshType<regionDynamicFvMesh>
{
    typedef dynamicFvMesh Type;
};


/*---------------------------------------------------------------------------*\
                  Class solverManagerBase::Control Declaration
\*---------------------------------------------------------------------------*/

typedef solutionControl Control;


/*---------------------------------------------------------------------------*\
                    Class solverManagerBase::Read Declaration
\*---------------------------------------------------------------------------*/

class Read
{
protected:

    // Protected Member Functions

        // Evolution

            //- Init sequence
            virtual void read() const = 0;


public:

    // Constructor

        Read() {};

    // Destructor

        virtual ~Read() {};
};


/*---------------------------------------------------------------------------*\
                   Class solverManagerBase::Init Declaration
\*---------------------------------------------------------------------------*/

class Init
{
private:

    // Private Data

        // Evolution

            //- Initalized?
            const bool initialized_;


protected:

    // Protected Member Functions

        // Evolution

            //- Init sequence
            virtual void init() const = 0;


public:

    // Constructor

        Init() : initialized_(false) {};

    // Destructor

        virtual ~Init() {};

    // Member Functions

        // Evolution

            //- Init wrapper
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
};


/*---------------------------------------------------------------------------*\
                 Class solverManagerBase::Settings Declaration
\*---------------------------------------------------------------------------*/

template <class MESH>
class Settings
:
    public virtual Read,
    public virtual Init
{
    //- Friends
    template<class M> friend class solverManager;
    template<class M> friend class Region;

private:

    // Private Member Functions

        // Evolution

            //- Init sequence (empty)
            virtual void init() const {};


public:

    // Constructor

        Settings() : Read(), Init() {};

    // Destructor

        virtual ~Settings() {};

    // Member Functions

        // Basic

            //- Mesh
            virtual MESH& mesh() const = 0;
};


/*---------------------------------------------------------------------------*\
                  Class solverManagerBase::Storage Declaration
\*---------------------------------------------------------------------------*/

template <class MESH>
class Storage
:
    public virtual Read,
    public virtual Init
{
    //- Friends
    template<class M> friend class solverManager;
    template<class M> friend class Region;

private:

    // Private Member Functions

        // Evolution

            //- Read sequence (empty)
            virtual void read() const {};


public:

    // Constructor

        Storage() : Read(), Init() {};

    // Destructor

        virtual ~Storage() {};

    // Member Functions

        // Basic

            //- Mesh
            virtual MESH& mesh() const = 0;
};


/*---------------------------------------------------------------------------*\
                  Class solverManagerBase::Region Declaration
\*---------------------------------------------------------------------------*/

template <class MESH>
class Region
:
    public virtual Read,
    public virtual Init
{
    //- Friends
    template<class M> friend class solverManager;
    template<class MM, class M> friend class Regions;

private:

    // Private Member Functions

        // Evolution

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
        virtual Settings<MESH>& settings() const = 0;

        //- Storage
        virtual Storage<MESH>& storage() const = 0;

        // Basic

            //- Mesh
            virtual MESH& mesh() const = 0;
};


/*---------------------------------------------------------------------------*\
                  Class solverManagerBase::Regions Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGERMESH, class MESH>
class Regions
:
    public virtual Read,
    public virtual Init
{
    //- Friends
    template<class M> friend class solverManager;

protected:

    // Protected typedefs

        typedef HashTable<Region<MESH>*> RegionHashTable;


private:

    // Private Data

        //- Hash table
        mutable autoPtr<RegionHashTable> regionHashTable_;

    // Private Member Functions

        // Evolution

            //- Read all regions in hash table
            virtual void read() const
            {
                forAllConstIter(
                    typename RegionHashTable, regionHashTable(), iter)
                {
                    iter()->read();
                }
            };
    
            //- Init all regions in hash table
            virtual void init() const
            {
                read();
    
                forAllConstIter(
                    typename RegionHashTable, regionHashTable(), iter)
                {
                    iter()->checkInit();
                }
            };


protected:

        // Protected Member Functions

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
// TODO: Remove?
        //- Default region
        virtual Region<MESH>& defaultRegion() const = 0;

    // Member Operators

        //- Find and return existing named region
        inline Region<MESH>& operator[](const word& name)
        {
            return *regionHashTable()[name];
        };

        // Basic

            //- Mesh
            virtual MANAGERMESH& mesh() const = 0;
};


/*---------------------------------------------------------------------------*\
                 Class solverManagerBase::Massages Declaration
\*---------------------------------------------------------------------------*/

//- Messages
template <class MESH>
class Messages
{
private:

    // Private Data

        //- Reference to args
        const argList& args_;

        //- Reference to time
        const Time& time_;

        //- Reference to mesh
        const MESH& mesh_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        Messages(const Messages&);

        //- Disallow default bitwise assignment
        void operator=(const Messages&);


public:

    // Constructor

        Messages
        (
            const argList& args,
            const Time& time,
            const MESH& mesh
        )
        :
            args_(args), time_(time), mesh_(mesh)
        {};

    // Destructor

        ~Messages() {};

    // Member Functions

        //- Info about starting time loop
        inline void newLine() const
        {
            Info << nl;
        }

        //- Info about starting time loop
        inline void startTimeLoop() const
        {
            Info << "Starting time loop"
                << endl;
        }

        //- Info about current time step
        inline void timeStep() const
        {
            Info << "Time step = "
                << time_.timeIndex()
                << endl;
        };

        //- Print current time step width
        inline void timeDeltaT() const
        {
            Info << "deltaT = "
                <<  time_.deltaT().value()
                << endl;
        };

        //- Print current time
        inline void timeIs() const
        {
            Info << "Time = "
                << time_.value()
                << endl;
        };

        //- Print execution time
        inline void executionTime() const
        {
            Info << "ExecutionTime = "
                << scalar(time_.elapsedCpuTime()) << " s"
                << endl;
        };

        //- Print end statement
        inline void end() const
        {
            Info << "End" << nl << endl;
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solverManagerBase

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
