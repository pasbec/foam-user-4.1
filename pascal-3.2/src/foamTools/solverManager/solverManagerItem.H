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
    solverManagerItem

Description
    ...

\*---------------------------------------------------------------------------*/

#ifndef solverManagerItem_H
#define solverManagerItem_H

#include "solverManagerWrapper.H"
#include "solverManagerStorage.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_ITEM_ADD_DECLARATION_START(WRAPPER, TYPE, NAME)                   \
                                                                              \
    public: class Item_##NAME                                                 \
    :                                                                         \
        public solverManagerWrappers::WRAPPER<ManagerType, TYPE>::Type        \
    {                                                                         \
        public: typedef typename                                              \
            solverManagerWrappers::WRAPPER<ManagerType, TYPE>::Type Wrapper;  \
                                                                              \
        public: typedef Item_##NAME ItemType;


# define SM_ITEM_ADD_DECLARATION_INTERFACES(NAME)                             \
                                                                              \
        private: const ManagerType& manager_;                                 \
        public: virtual const ManagerType& manager() const                    \
        {                                                                     \
            return manager_;                                                  \
        };                                                                    \
                                                                              \
        private: virtual void create(const word& init = "default") const;
        
                                                                              
# define SM_ITEM_ADD_DECLARATION_INTERFACES_REGION(NAME)                      \
                                                                              \
        public: virtual RegionType& region() const                            \
        {                                                                     \
            return static_cast<RegionType&>                                   \
            (                                                                 \
                manager().regions()[this->mesh().name()]                      \
            );                                                                \
        };                                                                    \
                                                                              \
        public: virtual RegionType::Control& control() const                  \
        {                                                                     \
            return region().control();                                        \
        };                                                                    \
                                                                              \
        public: virtual RegionType::Settings& settings() const                \
        {                                                                     \
            return region().settings();                                       \
        };                                                                    \
                                                                              \
        public: virtual StorageType& storage() const                          \
        {                                                                     \
            return region().storage();                                        \
        };


# define SM_ITEM_ADD_DECLARATION_INTERFACES_GLOBAL(NAME)                      \
                                                                              \
        public: virtual ManagerType::Control& control() const                 \
        {                                                                     \
            return manager().control();                                       \
        };                                                                    \
                                                                              \
        public: virtual ManagerType::Settings& settings() const               \
        {                                                                     \
            return manager().settings();                                      \
        };                                                                    \
                                                                              \
        public: virtual StorageType& storage() const                          \
        {                                                                     \
            return manager().storage();                                       \
        };


# define SM_ITEM_ADD_DECLARATION_CONSTRUCTION(NAME)                           \
                                                                              \
        private: Item_##NAME(const Item_##NAME&);                             \
        private: void operator=(const Item_##NAME&);                          \
                                                                              \
        public: Item_##NAME                                                   \
        (                                                                     \
            const ManagerType& manager,                                       \
            Wrapper::Mesh& mesh,                                              \
            const word& name                                                  \
        )                                                                     \
        :                                                                     \
            Wrapper(mesh, name),                                              \
            manager_(manager)                                                 \
        {};                                                                   \
                                                                              \
        public: virtual ~Item_##NAME() {};
        

# define SM_ITEM_ADD_DECLARATION_END()                                        \
                                                                              \
    };


# define SM_ITEM_ADD_DEFINITION(NAME)                                         \
                                                                              \
    private: mutable autoPtr<Item_##NAME> item_##NAME##Ptr_;                  \
                                                                              \
    public: Item_##NAME& item_##NAME() const                                  \
    {                                                                         \
        if (item_##NAME##Ptr_.empty())                                        \
        {                                                                     \
            SM_WRAPPER_MEMBER_SETINIT_P3(                                     \
                Item_##NAME, item_##NAME,                                     \
                this->manager(), this->mesh(), word(#NAME)                    \
            )                                                                 \
                                                                              \
            this->Base::itemHashTable().insert                                \
            (                                                                 \
                item_##NAME##Ptr_->name(),                                    \
                item_##NAME##Ptr_.operator->()                                \
            );                                                                \
        }                                                                     \
                                                                              \
        return item_##NAME##Ptr_();                                           \
    };


# define SM_ITEM_ADD_DATAACCESS(TYPE, NAME)                                   \
                                                                              \
    public: TYPE& NAME() const                                                \
    {                                                                         \
        return item_##NAME().ref();                                           \
    };                                                                        \
                                                                              \
    public: TYPE* NAME##Ptr() const                                           \
    {                                                                         \
        return item_##NAME().ptr();                                           \
    };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGIONITEM_ADD(TYPE, NAME)                                        \
                                                                              \
    SM_ITEM_ADD_DECLARATION_START(RegionItem, TYPE, NAME)                     \
    SM_ITEM_ADD_DECLARATION_INTERFACES()                                      \
    SM_ITEM_ADD_DECLARATION_INTERFACES_REGION()                               \
    SM_ITEM_ADD_DECLARATION_CONSTRUCTION(NAME)                                \
    SM_ITEM_ADD_DECLARATION_END()                                             \
    SM_ITEM_ADD_DEFINITION(NAME)                                              \
    SM_ITEM_ADD_DATAACCESS(TYPE, NAME)


# define SM_GLOBALITEM_ADD(TYPE, NAME)                                        \
                                                                              \
    SM_ITEM_ADD_DECLARATION_START(GlobalItem, TYPE, NAME)                     \
    SM_ITEM_ADD_DECLARATION_INTERFACES()                                      \
    SM_ITEM_ADD_DECLARATION_INTERFACES_GLOBAL()                               \
    SM_ITEM_ADD_DECLARATION_CONSTRUCTION(NAME)                                \
    SM_ITEM_ADD_DECLARATION_END()                                             \
    SM_ITEM_ADD_DEFINITION(NAME)                                              \
    SM_ITEM_ADD_DATAACCESS(TYPE, NAME)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solverManagerItem
{

/*---------------------------------------------------------------------------*\
                  Class solverManagerItem::Item Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER, class MESH, class TYPE>
class Item
:
    public solverManagerBase::Item<MESH, TYPE>,
    public solverManagerWrapper<MANAGER, MESH>
{
public:

    // Public typedefs

        typedef MANAGER Manager;

        typedef typename Manager::ManagerMesh ManagerMesh;
        typedef MESH Mesh;

        typedef solverManagerBase::Item<Mesh, TYPE> Base;
        typedef solverManagerWrapper<Manager, Mesh> Wrapper;

        typedef typename solverManagerBase::Control Control;
        typedef typename solverManagerBase::Settings<Mesh> Settings;
        typedef typename solverManagerBase::Storage<Mesh> Storage;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        Item(const Item&);

        //- Disallow default bitwise assignment
        void operator=(const Item&);

        // Basic

            //- Settings dictionary
            virtual const dictionary& settingsDict() const = 0;

            //- Storage dictionary
            virtual const dictionary& storageDict() const = 0;

        // Evolution

            //- Create sequence
            virtual void create(const word& init) const = 0;


public:

    // Constructor

        Item(Mesh& mesh, const word& name)
        :
            Base(name),
            Wrapper(mesh)
        {};

    // Destructor

        virtual ~Item() {};

    // Member Functions

        //- Manager
        virtual const Manager& manager() const = 0;

        //- Control
        virtual Control& control() const = 0;

        //- Settings
        virtual Settings& settings() const = 0;

        //- Storage
        virtual Storage& storage() const = 0;

        // Basic

            //- Mesh
            virtual Mesh& mesh() const
            {
                return this->Wrapper::wrappedMesh();
            }

        // Dictionaries

            //- Dictionary
            virtual const dictionary& dict() const
            {
                return storageDict();
            };
};


/*---------------------------------------------------------------------------*\
               Class solverManagerItem::RegionItem Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER, class MESH, class TYPE>
class RegionItem
:
    public Item<MANAGER, MESH, TYPE>
{
public:

    // Public typedefs

        typedef MANAGER Manager;

        typedef typename Manager::ManagerMesh ManagerMesh;
        typedef MESH Mesh;

        typedef Item<MANAGER, MESH, TYPE> Parent;
        typedef typename Parent::Base Base;
        typedef typename Parent::Wrapper Wrapper;

        typedef typename solverManagerBase::Region<Mesh> Region;

        typedef typename Parent::Control Control;
        typedef typename Parent::Settings Settings;
        typedef typename Parent::Storage Storage;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        RegionItem(const RegionItem&);

        //- Disallow default bitwise assignment
        void operator=(const RegionItem&);

        // Basic

            //- Settings dictionary
            virtual const dictionary& settingsDict() const
            {
                return this->Wrapper::settingsDict().subDict
                (
                    this->Parent::mesh().name()
                );
            };

            //- Storage dictionary
            virtual const dictionary& storageDict() const
            {
                return this->Wrapper::storageDict().subDict
                (
                    this->Parent::mesh().name()
                );
            };

        // Evolution

            //- Create sequence
            virtual void create(const word& init) const = 0;


public:

    // Constructor

        RegionItem(Mesh& mesh, const word& name)
        :
            Parent(mesh, name)
        {};

    // Destructor

        virtual ~RegionItem() {};

    // Member Functions

        //- Manager
        virtual const Manager& manager() const = 0;

        //- Region
        virtual Region& region() const = 0;

        //- Control
        virtual Control& control() const = 0;

        //- Settings
        virtual Settings& settings() const = 0;

        //- Storage
        virtual Storage& storage() const = 0;

};


/*---------------------------------------------------------------------------*\
              Class solverManagerItem::GlobalItem Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER, class MESH, class TYPE>
class GlobalItem
:
    public Item<MANAGER, MESH, TYPE>
{
public:

    // Public typedefs

        typedef MANAGER Manager;

        typedef typename Manager::ManagerMesh ManagerMesh;
        typedef MESH Mesh;

        typedef Item<MANAGER, MESH, TYPE> Parent;
        typedef typename Parent::Base Base;
        typedef typename Parent::Wrapper Wrapper;

        typedef typename Parent::Control Control;
        typedef typename Parent::Settings Settings;
        typedef typename Parent::Storage Storage;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        GlobalItem(const GlobalItem&);

        //- Disallow default bitwise assignment
        void operator=(const GlobalItem&);

        // Basic

            //- Settings dictionary
            virtual const dictionary& settingsDict() const
            {
                return this->Wrapper::settingsDict();
            };

            //- Storage dictionary
            virtual const dictionary& storageDict() const
            {
                return this->Wrapper::storageDict();
            };

        // Evolution

            //- Create sequence
            virtual void create(const word& init) const = 0;


public:

    // Constructor

        GlobalItem(Mesh& mesh, const word& name)
        :
            Parent(mesh, name)
        {};

    // Destructor

        virtual ~GlobalItem() {};

    // Member Functions

        //- Manager
        virtual const Manager& manager() const = 0;

        //- Control
        virtual Control& control() const = 0;

        //- Settings
        virtual Settings& settings() const = 0;

        //- Storage
        virtual Storage& storage() const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solverManagerItem

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solverManagerWrappers
{
    template <class MANAGER, class TYPE>
    struct RegionItem
    {
        typedef MANAGER Manager;

        typedef TYPE ItemType;

        typedef typename MANAGER::ManagerMesh ManagerMesh;
        typedef typename MANAGER::Mesh Mesh;

        typedef solverManagerItem::RegionItem<Manager, Mesh, ItemType> Type;
    };


    template <class MANAGER, class TYPE>
    struct GlobalItem
    {
        typedef MANAGER Manager;

        typedef TYPE ItemType;

        typedef typename MANAGER::ManagerMesh ManagerMesh;
        typedef typename MANAGER::ManagerMesh Mesh;

        typedef solverManagerItem::GlobalItem<Manager, Mesh, ItemType> Type;
    };
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
