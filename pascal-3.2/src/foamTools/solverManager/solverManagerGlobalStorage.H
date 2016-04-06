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
    solverManagerGlobalStorage

Description
    ...

\*---------------------------------------------------------------------------*/

#ifndef solverManagerGlobalStorage_H
#define solverManagerGlobalStorage_H

#include "solverManagerWrapper.H"
#include "solverManagerStorage.H"

// TODO: makeOptional_ has to be connected to settings!

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_GLOBALSTORAGE_NEW(MANAGER)                                        \
                                                                              \
    SM_WRAPPER_NEW(GlobalStorage, Storage, MANAGER, Storage)                 \
                                                                              \
    private: virtual void init() const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_GLOBALSTORAGE_MEMBER_ADD(TYPE, NAME)                              \
                                                                              \
        SM_STORAGE_MEMBER_ADD(TYPE, NAME)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_GLOBALSTORAGE_ADD()                                               \
                                                                              \
    SM_WRAPPER_END()                                                          \
                                                                              \
    SM_ADD_GLOBALSTORAGE()


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class solverManagerGlobalStorage Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER>
class solverManagerGlobalStorage
:
    public solverManagerCore<typename MANAGER::ManagerMesh>::Storage,
    public solverManagerWrapper<MANAGER, typename MANAGER::ManagerMesh>
{
public:

    // Public typedefs

        typedef MANAGER Manager;

        typedef typename Manager::ManagerMesh ManagerMesh;
        typedef typename Manager::ManagerMesh Mesh;

        typedef solverManagerWrapper<Manager, Mesh> Wrapper;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        solverManagerGlobalStorage(const solverManagerGlobalStorage&);

        //- Disallow default bitwise assignment
        void operator=(const solverManagerGlobalStorage&);

public:

    // Constructor

        solverManagerGlobalStorage(const Manager& manager, Mesh& mesh)
        :
            Wrapper(manager, mesh)
        {};

    // Destructor

        virtual ~solverManagerGlobalStorage() {};

    // Member Functions

        // Basic

            //- Mesh
            inline Mesh& mesh() const
            {
                return this->Wrapper::wrappedMesh();
            }

        // Dictionaries

            //- Dictionary
            inline const dictionary& dict() const
            {
                return this->Wrapper::storageDict();
            };


            //- Read dictionary of named item in storage (if present)
            bool itemDict
            (
                const word& name,
                dictionary& itemDict
            ) const
            {
                if
                (
                    dict().found(name)
                 && dict().isDict(name)
                )
                {
                    itemDict = dict().subDict(name);

                    Switch enabled;

                    if (itemDict.readIfPresent("enabled", enabled))
                    {
                        return enabled;
                    }
                }

                return false;
            };
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
