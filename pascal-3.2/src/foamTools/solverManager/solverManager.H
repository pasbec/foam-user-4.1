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

\*---------------------------------------------------------------------------*/

#ifndef solverManager_H
#define solverManager_H


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGIONSCOPE_MESHONLY(MANAGER, SCOPETYPE, SCOPE)                   \
                                                                              \
    MANAGER::SCOPETYPE::Mesh& mesh =                                          \
        manager.regions().SCOPE().mesh();                                     \
                                                                              \
    if (&mesh != &mesh) {};


# define SM_REGIONSCOPE_NOMESH(MANAGER, SCOPETYPE, SCOPE)                     \
                                                                              \
    MANAGER::ManagerMesh& globalMesh = manager.mesh();                        \
    MANAGER::Control& globalControl = manager.control();                      \
    MANAGER::Settings& globalSettings = manager.settings();                   \
    MANAGER::Storage& globalStorage = manager.storage();                      \
                                                                              \
    if (&globalMesh != &globalMesh) {};                                       \
    if (&globalControl != &globalControl) {};                                 \
    if (&globalSettings != &globalSettings) {};                               \
    if (&globalStorage != &globalStorage) {};                                 \
                                                                              \
    MANAGER::SCOPETYPE::Control& control =                                    \
        manager.regions().SCOPE().control();                                  \
    MANAGER::SCOPETYPE::Settings& settings =                                  \
        manager.regions().SCOPE().settings();                                 \
    MANAGER::SCOPETYPE::Storage& storage =                                    \
        manager.regions().SCOPE().storage();                                  \
                                                                              \
    if (&control != &control) {};                                             \
    if (&settings != &settings) {};                                           \
    if (&storage != &storage) {};


# define SM_REGIONSCOPE(MANAGER, SCOPETYPE, SCOPE)                            \
                                                                              \
    SM_REGIONSCOPE_NOMESH(MANAGER, SCOPETYPE, SCOPE)                          \
    SM_REGIONSCOPE_MESHONLY(MANAGER, SCOPETYPE, SCOPE)


# define SM_REGIONSCOPE_DEFAULTREGION(MANAGER)                                \
                                                                              \
    SM_REGIONSCOPE_NOMESH(MANAGER, DefaultRegion, defaultRegion);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solverManagerCore.H"
#include "solverManagerWrapper.H"
#include "solverManagerGlobalSettings.H"
#include "solverManagerGlobalStorage.H"
#include "solverManagerSettings.H"
#include "solverManagerStorage.H"
#include "solverManagerRegion.H"
#include "solverManagerRegions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
