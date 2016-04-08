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
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template <class MESH> class solverManager;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solverManagerBase
{

template <class MESH> class Mesh;
template <class MESH> class Settings;
template <class MESH, class TYPE> class Item;
template <class MESH> class Storage;
template <class MESH> class Region;
template <class MANAGERMESH, class MESH> class Regions;

/*---------------------------------------------------------------------------*\
                  Class solverManagerBase::Mesh Declaration
\*---------------------------------------------------------------------------*/

template <typename MESH>
struct Mesh
{
    typedef typename MESH::Mesh Type;
};

template <>
struct Mesh<dynamicFvMesh>
{
    typedef dynamicFvMesh Type;
};

template <>
struct Mesh<regionDynamicFvMesh>
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
private:

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
            mutable Switch initialized_;

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

            //- Init sequence wrapper
            bool checkInit() const
            {
                if (!initialized_)
                {
                    initialized_ = true;

                    init();
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
    public Read,
    public Init
{
    //- Friends
    template<class M> friend class solverManager;
    template<class M> friend class Region;

private:

    // Private Member Functions

        // Evolution

            //- Read sequence
            virtual void read() const = 0;

            //- Init sequence (empty)
            virtual void init() const {};


public:

    // Constructor

        Settings() : Read(), Init() {};

    // Destructor

        virtual ~Settings() {};

    // Member Functions

        //- Control
        virtual Control& control() const = 0;

        //- Storage
        virtual Storage<MESH>& storage() const = 0;

        // Basic

            //- Mesh
            virtual MESH& mesh() const = 0;
};


/*---------------------------------------------------------------------------*\
                   Class solverManagerBase::Item Declaration
\*---------------------------------------------------------------------------*/

class ItemCore
:
    public Read,
    public Init
{
    //- Friends
    template<class M, class T> friend class Item;
    template<class M> friend class Storage;

private:

    // Private Data

        //- Name
        word name_;

        //- Disabled switch
        mutable Switch disabled_;

        //- Locked switch
        mutable Switch locked_;

        //- Optional switch
        mutable Switch optional_;

    // Private Member Functions

        // Evolution

            //- Read sequence (empty)
            virtual void read() const {};

            //- Init sequence
            virtual void init() const
            {
                if (!disabled())
                {
                    checkMake();
                }
                else
                {
                    FatalErrorIn
                    (
                        "solverManagerBase::ItemCore::"
                        + this->ItemCore::name() + "::init()"
                    )
                        << "Item is globally disabled and cannot be used."
                            << " Check your settings!"
                            << abort(FatalError);

                }
            };

            //- Make sequence wrapper
            virtual void checkMake() const = 0;


protected:

    // Protected Member Functions

        // Disabling

            //- Reference to disabled switch
            Switch& disabledSwitch() const
            {
                return disabled_;
            };

        // Locking

            //- Reference to locked switch
            Switch& lockedSwitch() const
            {
                return locked_;
            };

            //- Is item locked?
            inline bool locked() const
            {
                return locked_;
            };

            //- Lock item
            inline void lock() const
            {
                locked_ = true;
            };

            //- Unlock item
            inline void unlock() const
            {
                locked_ = false;
            };


public:

    // Constructor

        ItemCore(const word& name)
        :
            Read(),
            Init(),
            name_(name),
            disabled_(false),
            locked_(false),
            optional_(false)
        {};

    // Destructor

        virtual ~ItemCore() {};

    // Member Functions

        //- Name
        inline const word& name() const
        {
            return name_;
        };

        // Pointer

            //- Item is empty
            virtual bool empty() const = 0;

            //- Item is valid
            virtual bool valid() const = 0;

        // Disabling

            //- Is item disabled?
            inline bool disabled() const
            {
                return disabled_;
            };

            //- Disable item
            virtual void disable() const = 0;

            //- Enable item
            virtual void enable() const = 0;

        // Optional

            //- Is item optional?
            inline bool optional() const
            {
                return optional_;
            };

            //- Make item optional
            inline void makeOptional() const
            {
                optional_ = true;
            };

            //- Make item mandatory
            inline void makeMandatory() const
            {
                optional_ = false;
            };
};


template <class MESH, class TYPE>
class Item
:
    public ItemCore
{
private:

    // Private Data

        //- Item data
        mutable autoPtr<TYPE> itemPtr_;

    // Private Member functions

        // Evolution

            //- Create sequence
            virtual void create(const word& init) const = 0;

            //- Make sequence
            void make() const
            {
                if (!this->ItemCore::locked())
                {
                    this->ItemCore::lockedSwitch() = true;

                    Info << "Create " << this->ItemCore::name()
                        << " for region " << mesh().name()
                        << endl;

                    word init = dict().lookupOrDefault(name(), word("default"));

                    create(init);

                    this->ItemCore::lockedSwitch() = false;
                }
                else
                {
                    FatalErrorIn
                    (
                        "solverManagerBase::Item::"
                        + this->ItemCore::name() + "::make()"
                    )
                        << "Pointer for item is part of a circular "
                        << "reference. Check all occurances!"
                        << abort(FatalError);
                }
            };

            //- Make sequence wrapper
            virtual void checkMake() const
            {
                if (itemPtr_.empty())
                {
                    make();

                    if (itemPtr_.empty())
                    {
                        FatalErrorIn
                        (
                            "solverManagerBase::Item::"
                            + this->ItemCore::name() + "::checkMake()"
                        )
                            << "Pointer for item is still unallocated "
                            << "after creation. Make sure you have set "
                            << "the pointer via set(TYPE*) in its "
                            << "create-function!"
                            << abort(FatalError);
                    }
                }
            };


protected:

    // Private Member functions

        // Access

            //- Return auto pointer
            inline autoPtr<TYPE> get() const
            {
                return itemPtr_;
            };

            //- Set auto pointer
            inline void set(TYPE* p) const
            {
                itemPtr_.set(p);
            };


public:

    // Constructor

        Item(const word& name) : ItemCore(name) {};

    // Destructor

        virtual ~Item() {};

    // Member Functions

        //- Control
        virtual Control& control() const = 0;

        //- Settings
        virtual Settings<MESH>& settings() const = 0;

        //- Storage
        virtual Storage<MESH>& storage() const = 0;

        // Basic

            //- Mesh
            virtual MESH& mesh() const = 0;

        // Dictionaries

            virtual const dictionary& dict() const = 0;

        // Access

            //- Return auto pointer (DANGEROUS!)
            inline autoPtr<TYPE> aPtr() const
            {
                return itemPtr_;
            };

            //- Return raw pointer (DANGEROUS!)
            inline TYPE* rPtr() const
            {
                return itemPtr_.operator->();
            };

            //- Return pointer
            inline TYPE* ptr() const
            {
                this->ItemCore::init();

                return itemPtr_.operator->();
            };

            //- Return reference
            inline TYPE& ref() const
            {
                this->ItemCore::init();

                return itemPtr_();
            };

        // Pointer

            //- Item is empty
            virtual bool empty() const
            {
                return itemPtr_.empty();
            };

            //- Item is valid
            virtual bool valid() const
            {
                return itemPtr_.valid();
            };

        // Disabling

            //- Disable item
            virtual void disable() const
            {
                if (!this->ItemCore::disabled())
                {
                    if (valid())
                    {
                        FatalErrorIn
                        (
                            "solverManagerBase::Item::"
                            + this->ItemCore::name() + "::disable()"
                        )
                            << "Item is still allocated and cannot "
                            <<"be disabled. Clear it first!"
                            << abort(FatalError);
                    }
                    else
                    {
                        this->ItemCore::disabledSwitch() = true;
                    }
                }
            };

            //- Enable item
            virtual void enable() const
            {
                if (this->ItemCore::disabled())
                {
                    FatalErrorIn
                    (
                        "solverManagerBase::Item::"
                        + this->ItemCore::name() + "::enable()"
                    )
                        << "Item is already enabled!"
                        << abort(FatalError);
                }
                else
                {
                    this->ItemCore::disabledSwitch() = false;
                }
            };
};


/*---------------------------------------------------------------------------*\
                  Class solverManagerBase::Storage Declaration
\*---------------------------------------------------------------------------*/

template <class MESH>
class Storage
:
    public Read,
    public Init
{
    //- Friends
    template<class M> friend class solverManager;
    template<class M> friend class Region;

protected:

    // Protected typedefs

        typedef HashTable<ItemCore*> ItemHashTable;


private:

    // Private Data

        //- Hash table
        mutable autoPtr<ItemHashTable> itemHashTable_;

    // Private Member Functions

        // Evolution

            //- Read sequence (empty)
            virtual void read() const {};

            //- Init all regions in hash table
            virtual void init() const
            {
                read();

                create();

                forAllConstIter(
                    typename ItemHashTable, itemHashTable(), iter)
                {
                    iter()->checkInit();
                }
            };

            //- Create sequence
            virtual void create() const = 0;


protected:

        // Protected Member Functions

            //- Hash table
            inline ItemHashTable& itemHashTable() const
            {
                if (itemHashTable_.empty())
                {
                    itemHashTable_.set(new ItemHashTable());
                }

                return itemHashTable_();
            };

public:

    // Constructor

        Storage() : Read(), Init() {};

    // Destructor

        virtual ~Storage() {};

    // Member Functions

        //- Control
        virtual Control& control() const = 0;

        //- Settings
        virtual Settings<MESH>& settings() const = 0;

        //- Find existing named item and return ref cast to T
        template<class T>
        inline Item<MESH, T>& lookup(const word& name)
        {
            return *dynamic_cast<Item<MESH, T>* >(itemHashTable()[name]);
        };

        // Basic

            //- Mesh
            virtual MESH& mesh() const = 0;

    // Member Operators

        //- Find and return existing named item
        inline ItemCore& operator[](const word& name)
        {
            return *itemHashTable()[name];
        };
};


/*---------------------------------------------------------------------------*\
                  Class solverManagerBase::Region Declaration
\*---------------------------------------------------------------------------*/

template <class MESH>
class Region
:
    public Read,
    public Init
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

        //- Control
        virtual Control& control() const = 0;

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
    public Read,
    public Init
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

        //- Find existing named region and ref cast to T
        template<class T>
        inline T& lookup(const word& name)
        {
            return *dynamic_cast<T*>(regionHashTable()[name]);
        };

        // Basic

            //- Mesh
            virtual MANAGERMESH& mesh() const = 0;

    // Member Operators

        //- Find and return existing named region
        inline Region<MESH>& operator[](const word& name)
        {
            return *regionHashTable()[name];
        };
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
