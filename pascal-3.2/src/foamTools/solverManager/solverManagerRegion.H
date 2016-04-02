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
    solverManagerRegion

Description
    ...

\*---------------------------------------------------------------------------*/

#ifndef solverManagerRegion_H
#define solverManagerRegion_H

#include "solverManager.H"
#include "solverManagerBaseNested.H"
#include "solverManagerStorage.H"

// TODO: Macro new for regionMeshes where mesh is the mesh of one region:
//       Here and in storage MeshType should be used for all "mesh"-references
//       with a template specification for dydnamicFvMesh

// TODO: Use mesh().name() for the virtual functions settingsRegionDict(...), storageRegionDict(...)
//       and storageItemDict(...) instead of polyMesh::defaultRegion as default
//       parameter for the region word

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define smRegion_new(MANAGER, NAME)                                          \
                                                                              \
    smBaseNested_new(Region, solverManagerRegion, MANAGER, NAME)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define smRegion_end()                                                       \
                                                                              \
    smBaseNested_end()


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define smRegion_addControl_mesh(TYPE)                                       \
                                                                              \
    public: typedef TYPE Control;                                             \
                                                                              \
    smBaseNested_addModItemNoInit_P1(Control, control, virtual, this->mesh())


# define smRegion_addSettings_manager()                                       \
                                                                              \
    smBaseNested_addVirtualItem_P1(Settings, settings, this->manager())


# define smRegion_addStorage_manager()                                        \
                                                                              \
    smBaseNested_addVirtualItem_P1(Storage, storage, this->manager())


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class solverManagerRegion Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER>
class solverManagerRegion
:
    public solverManagerBaseNested<MANAGER>
{
public:

    // Public typedefs

        typedef MANAGER Manager;
        typedef typename MANAGER::ManagerMesh ManagerMesh;
        typedef typename ManagerMesh::Mesh ThisMesh;

        typedef solutionControl Control;
        typedef solverManagerStorage<MANAGER> Storage;


private:

//     // Private Data
//
//         //- Reference to mesh
//         ThisMesh& mesh_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        solverManagerRegion(const solverManagerRegion&);

        //- Disallow default bitwise assignment
        void operator=(const solverManagerRegion&);


protected:

//     // Protected Member Functions
//
//         // Basic
//
//             //- Return reference to mesh
//             ThisMesh& mesh() const
//             {
//                 return mesh_;
//             }

//             //- Reference to settings dictionary for named region
//             const dictionary& settingsRegionDict
//             (
//                 const word& region = ????????
//             ) const
//             {
//                 return manager_.settingsRegionDict(region);
//             };
//
//             //- Reference to storage dictionary for named region
//             const dictionary& storageRegionDict
//             (
//                 const word& region = ????????
//             ) const
//             {
//                 return manager_.storageRegionDict(region);
//             };
//
//             //- Reference to dictionary of named item of storage
//             bool storageItemDict
//             (
//                 const word& item,
//                 dictionary& dict,
//                 const word& region = ????????
//             ) const
//             {
//                 return manager_.storageItemDict(item, dict, region);
//             };


public:

    // Constructor

        solverManagerRegion(const MANAGER& manager)
        :
            solverManagerBaseNested<MANAGER>(manager)
        {};

    // Destructor

        virtual ~solverManagerRegion() {};
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
