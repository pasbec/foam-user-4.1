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

#ifndef solverManagerItem__H
#define solverManagerItem__H

#include "solverManagerWrapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_ITEM_ADD_DECLARATION_START(WRAPPER, TYPE, NAME)                   \
                                                                              \
    public: class Item_##NAME                                                 \
    :                                                                         \
        public solverManagerWrappers::WRAPPER<Manager, TYPE>                  \
    {                                                                         \
        public: typedef typename                                              \
            solverManagerWrappers::WRAPPER<Manager, TYPE> Wrapper;            \
        public: typedef typename                                              \
            solverManagerWrappers::WRAPPER<Manager, TYPE>::Mesh Mesh;


# define SM_ITEM_ADD_DECLARATION_INTERFACES(NAME)                             \
                                                                              \
        private: const Manager& manager_;                                     \
        public: virtual const Manager& manager() const                        \
        {                                                                     \
            return manager_;                                                  \
        };                                                                    \
                                                                              \
        private: virtual void create(const word& ccase = "default") const;


# define SM_ITEM_ADD_DECLARATION_INTERFACES_REGION(NAME)                      \
                                                                              \
        public: inline Manager::Control& globalControl() const                \
        {                                                                     \
            return manager().control();                                       \
        };                                                                    \
                                                                              \
        public: inline Manager::Settings& globalSettings() const              \
        {                                                                     \
            return manager().settings();                                      \
        };                                                                    \
                                                                              \
        public: inline Manager::Storage& globalStorage() const                \
        {                                                                     \
            return manager().storage();                                       \
        };                                                                    \
                                                                              \
        public: virtual ThisRegion& region() const                            \
        {                                                                     \
            return static_cast<ThisRegion&>                                   \
            (                                                                 \
                manager().regions()[regionLabel]                              \
            );                                                                \
        };                                                                    \
                                                                              \
        public: inline ThisRegion::Control& control() const                   \
        {                                                                     \
            return region().control();                                        \
        };                                                                    \
                                                                              \
        public: inline ThisRegion::Settings& settings() const                 \
        {                                                                     \
            return region().settings();                                       \
        };                                                                    \
                                                                              \
        public: virtual ThisRegion::Storage& storage() const                  \
        {                                                                     \
            return region().storage();                                        \
        };


# define SM_ITEM_ADD_DECLARATION_INTERFACES_GLOBAL(NAME)                      \
                                                                              \
        public: inline Manager::Control& control() const                      \
        {                                                                     \
            return manager().control();                                       \
        };                                                                    \
                                                                              \
        public: inline Manager::Settings& settings() const                    \
        {                                                                     \
            return manager().settings();                                      \
        };                                                                    \
                                                                              \
        public: virtual Manager::Storage& storage() const                     \
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
            const Manager& manager,                                           \
            Mesh& mesh,                                                       \
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
    SM_WRAPPER_MEMBER_ADD_P3(                                                 \
        Item_##NAME, item_##NAME, virtual,                                    \
        this->manager(), this->mesh(), word(#NAME))


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
    SM_ITEM_ADD_DECLARATION_START(RegionItem_, TYPE, NAME)                    \
    SM_ITEM_ADD_DECLARATION_INTERFACES()                                      \
    SM_ITEM_ADD_DECLARATION_INTERFACES_REGION()                               \
    SM_ITEM_ADD_DECLARATION_CONSTRUCTION(NAME)                                \
    SM_ITEM_ADD_DECLARATION_END()                                             \
    SM_ITEM_ADD_DEFINITION(NAME)                                              \
    SM_ITEM_ADD_DATAACCESS(TYPE, NAME)


# define SM_GLOBALITEM_ADD(TYPE, NAME)                                        \
                                                                              \
    SM_ITEM_ADD_DECLARATION_START(GlobalItem_, TYPE, NAME)                    \
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

namespace solverManagerWrappers
{

/*---------------------------------------------------------------------------*\
                  Class solverManagerWrappers::Item_ Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER, class MESH, class TYPE>
class Item_
:
    public solverManagerBase::Item_<MESH, TYPE>,
    public solverManagerWrappers::Base<MANAGER, MESH>
{
public:

    // Public typedefs

        typedef MESH Mesh;


private:

    // Private typedefs

        typedef solverManagerBase::Item_<MESH, TYPE> Base;
        typedef solverManagerWrappers::Base<MANAGER, MESH> Wrapper;

        typedef typename solverManagerBase::Storage<MESH> Storage;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        Item_(const Item_&);

        //- Disallow default bitwise assignment
        void operator=(const Item_&);

        // Basic

            //- Settings dictionary
            virtual const dictionary& settingsDict() const = 0;

            //- Storage dictionary
            virtual const dictionary& storageDict() const = 0;

        // Evolution

            //- Create sequence
            virtual void create(const word& ccase) const = 0;


public:

    // Constructor

        Item_(MESH& mesh, const word& name)
        :
            Base(name),
            Wrapper(mesh)
        {};

    // Destructor

        virtual ~Item_() {};

    // Member Functions

        //- Storage
        virtual Storage& storage() const = 0;

        //- Manager
        virtual const MANAGER& manager() const = 0;

        // Basic

            //- Mesh
            virtual MESH& mesh() const
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
               Class solverManagerItem::RegionItem_ Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER, class TYPE>
class RegionItem_
:
    public Item_<MANAGER, typename MANAGER::Mesh, TYPE>
{
public:

    // Public typedefs

        typedef typename MANAGER::Mesh Mesh;


private:

    // Private typedefs

        typedef Item_<MANAGER, Mesh, TYPE> Parent;
        typedef solverManagerBase::Item_<Mesh, TYPE> Base;
        typedef solverManagerWrappers::Base<MANAGER, Mesh> Wrapper;

        typedef solverManagerBase::Region_
        <
            typename MANAGER::ManagerMesh,
            Mesh,
            MANAGER::nRegions
        > Region_;
        typedef typename solverManagerBase::Storage<Mesh> Storage;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        RegionItem_(const RegionItem_&);

        //- Disallow default bitwise assignment
        void operator=(const RegionItem_&);

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
            virtual void create(const word& ccase) const = 0;


public:

    // Constructor

        RegionItem_(Mesh& mesh, const word& name)
        :
            Parent(mesh, name)
        {};

    // Destructor

        virtual ~RegionItem_() {};

    // Member Functions

        //- Manager
        virtual const MANAGER& manager() const = 0;

        //- Region
        virtual Region_& region() const = 0;

        //- Storage
        virtual Storage& storage() const = 0;

};


/*---------------------------------------------------------------------------*\
              Class solverManagerItem::GlobalItem_ Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER, class TYPE>
class GlobalItem_
:
    public Item_<MANAGER, typename MANAGER::ManagerMesh, TYPE>
{
public:

    // Public typedefs

        typedef typename MANAGER::ManagerMesh Mesh;


private:

    // Private typedefs

        typedef Item_<MANAGER, Mesh, TYPE> Parent;
        typedef solverManagerBase::Item_<Mesh, TYPE> Base;
        typedef solverManagerWrappers::Base<MANAGER, Mesh> Wrapper;

        typedef typename solverManagerBase::Storage<Mesh> Storage;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        GlobalItem_(const GlobalItem_&);

        //- Disallow default bitwise assignment
        void operator=(const GlobalItem_&);

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
            virtual void create(const word& ccase) const = 0;


public:

    // Constructor

        GlobalItem_(Mesh& mesh, const word& name)
        :
            Parent(mesh, name)
        {};

    // Destructor

        virtual ~GlobalItem_() {};

    // Member Functions

        //- Manager
        virtual const MANAGER& manager() const = 0;

        //- Storage
        virtual Storage& storage() const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solverManagerWrappers

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
