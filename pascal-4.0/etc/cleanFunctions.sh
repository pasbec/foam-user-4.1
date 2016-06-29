#!/bin/bash

#---------------------------------*- sh -*-------------------------------------
# =========                 |
# \\      /  F ield         | foam-extend: Open Source CFD
#  \\    /   O peration     |
#   \\  /    A nd           | For copyright notice see file Copyright
#    \\/     M anipulation  |
#------------------------------------------------------------------------------
# License
#     This file is part of foam-extend.
#
#     foam-extend is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     foam-extend is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.
#
# Script
#     CleanFunctions
#
# Description
#
#------------------------------------------------------------------------------

cleanPolyMesh ()
{
    baseDir='constant'
    [[ ! -z "$1" ]] && baseDir="${baseDir}/$1"

    rm -rf "$baseDir/polyMesh/allOwner"* \
           "$baseDir/polyMesh/cell"* \
           "$baseDir/polyMesh/face"* \
           "$baseDir/polyMesh/meshModifiers"* \
           "$baseDir/polyMesh/owner"* \
           "$baseDir/polyMesh/neighbour"* \
           "$baseDir/polyMesh/point"* \
           "$baseDir/polyMesh/edge"* \
           "$baseDir/polyMesh/zoneToPatchName" \
           "$baseDir/cellLevel" \
           "$baseDir/polyMesh/cellLevel"* \
           "$baseDir/pointLevel" \
           "$baseDir/polyMesh/pointLevel"* \
           "$baseDir/polyMesh/refinementHistory"* \
           "$baseDir/polyMesh/surfaceIndex"* \
           "$baseDir/cellToRegion" \
           "$baseDir/polyMesh/sets/" \
           > /dev/null 2>&1
}

cleanFaMesh ()
{
    baseDir='constant'
    [[ ! -z "$1" ]] && baseDir="${baseDir}/$1"

    rm -rf "$baseDir/faMesh/faceLabels"* > /dev/null 2>&1
    rm -rf "$baseDir/faMesh/faBoundary"* > /dev/null 2>&1
}

cleanTimeDirectories ()
{
    [[ "$1" != 'silent' ]] && echo "Cleaning $PWD time directories"

    nZeros=0
    zeros=""
    while [ $nZeros -lt 8 ] ; do
        timeDir="0.${zeros}[1-9]*"
        rm -rf ${timeDir} > /dev/null 2>&1
        rm -rf ./-${timeDir} > /dev/null 2>&1
        zeros=`printf %0${nZeros}d 0`
        nZeros=$(($nZeros + 1))
    done
    rm -rf ./[1-9]* ./-[1-9]* > /dev/null 2>&1
}

cleanTimeDirectoriesParallel ()
{
    [[ "$1" != 'silent' ]] && echo "Cleaning $PWD time directories"

    nZeros=0
    zeros=""
    while [ $nZeros -lt 8 ] ; do
        timeDir="0.${zeros}[1-9]*"
        rm -rf processor*/${timeDir} > /dev/null 2>&1
        rm -rf processor*/-${timeDir} > /dev/null 2>&1
        zeros=`printf %0${nZeros}d 0`
        nZeros=$(($nZeros + 1))
    done
    rm -rf processor*/[1-9]* processor*/-[1-9]* > /dev/null 2>&1
}

cleanCase ()
{
    if [[ -z "$1" ]]; then
        echo "Cleaning case"
    else
        echo "Cleaning region $1"
    fi

    cleanPolyMesh "$1"
    cleanFaMesh "$1"

    if [[ -z "$1" ]]; then
        cleanTimeDirectories 'silent'

        rm -rf processor* > /dev/null 2>&1
        rm -rf probes* > /dev/null 2>&1
        rm -rf forces* > /dev/null 2>&1

        rm -rf system/machines \
            VTK \
            > /dev/null 2>&1

        rm -rf ./log \
            ./log.* \
            ./log-* \
            ./logSummary.* \
            ./.setSet \
            ./.fxLock \
            ./*.xml \
            ./ParaView* \
            ./paraFoam* \
            ./*.blockMesh \
            ./*.OpenFOAM \
            ./*.foam \
            ./*.sMesh \
            ./PlyParser_FoamFileParser_parsetab.py \
            > /dev/null 2>&1

#         for f in $(find . -name "*Dict")
#         do
#             sed -e /arguments/d $f > temp.$$
#             mv temp.$$ $f
#         done
    fi
}


#####


cleanToolsComsol ()
{
    rm -rf 'tools/comsol/log' > /dev/null 2>&1
    rm -rf 'tools/comsol/tmp' > /dev/null 2>&1
    rm -rf 'tools/comsol/'*'.class' > /dev/null 2>&1
    rm -rf 'tools/comsol/'*'.class.status' > /dev/null 2>&1
}

cleanTools ()
{
    cleanToolsComsol

    rm -rf 'tools/'*'.pyc' \
           'PlyParser_FoamFileParser_parsetab.py' \
            > /dev/null 2>&1
}


#####


wipePolyMesh ()
{
    baseDir='constant'
    [[ ! -z "$1" ]] && baseDir="${baseDir}/$1"

    rm -rf "$baseDir/polyMesh" > /dev/null 2>&1
}

wipeFaMesh ()
{
    baseDir='constant'
    [[ ! -z "$1" ]] && baseDir="${baseDir}/$1"

    rm -rf "$baseDir/faMesh" > /dev/null 2>&1
}

wipeZero ()
{
    rm -rf '0' > /dev/null 2>&1
}

wipeFaMeshParallel ()
{
    for p in $(ls -1d 'processor'*); do

        baseDir="$p/constant"
        [[ ! -z "$1" ]] && baseDir="${baseDir}/$1"

        rm -rf "$baseDir/faMesh" > /dev/null 2>&1

    done
}


#####


wipeOrgCopyPolyMesh ()
{
    baseDir='constant'
    [[ ! -z "$1" ]] && baseDir="${baseDir}/$1"

    [[ -d "$baseDir/polyMesh.org" ]] && wipePolyMesh "$1"
}

wipeOrgCopyFaMesh ()
{
    baseDir='constant'
    [[ ! -z "$1" ]] && baseDir="${baseDir}/$1"

    [[ -d "$baseDir/faMesh.org" ]] && wipeFaMesh "$1"
}

wipeOrgCopyZero ()
{
    [[ -d '0.org' ]] && wipeZero
}

wipeOrgCopy ()
{
    wipeOrgCopyPolyMesh "$1"
    wipeOrgCopyFaMesh "$1"

    [[ -z "$1" ]] && wipeOrgCopyZero
}

# DEPRECATED
wipeOrg ()
{
    echo "[DEPRECATED]: The 'wipeOrg'-function is deprecated as its name is misleading."
    echo "              Please use the new 'wipeOrgCopy'-function instead!"
    wipeOrgCopy "$1"
}


#####


restorePolyMeshOrg ()
{
    baseDir='constant'
    [[ ! -z "$1" ]] && baseDir="${baseDir}/$1"

    wipePolyMesh "$1"
    mkdir "$baseDir/polyMesh" > /dev/null 2>&1
    cp -r "$baseDir/polyMesh.org/"* 'constant/polyMesh' > /dev/null 2>&1
}

restoreZeroOrg ()
{
    wipeZero
    mkdir '0' > /dev/null 2>&1
    cp -r "0.org/"* '0' > /dev/null 2>&1
}

restoreOrg ()
{
    restorePolyMeshOrg "$1"
    restoreFaMeshOrg "$1"

    [[ -z "$1" ]] && restoreZeroOrg
}


####


restoreFaMeshOrg ()
{
    wipeFaMesh "$1"

    baseDir='constant'
    baseDirConstant='constant'
    if [[ ! -z "$1" ]]; then
        baseDir="${baseDir}/$1"
        baseDirConstant="${baseDirConstant}/$1"
    fi



    [[ ! -d "$baseDir/faMesh" ]] && mkdir "$baseDir/faMesh" > /dev/null 2>&1
    cp -r "$baseDirConstant/faMesh.org/"* "$baseDir/faMesh" > /dev/null 2>&1
}

restoreFaMeshOrgParallel ()
{
    wipeFaMeshParallel "$1"

    for p in $(ls -1d 'processor'*); do

        baseDir="$p/constant"
        baseDirConstant='constant'
        if [[ ! -z "$1" ]]; then
            baseDir="${baseDir}/$1"
            baseDirConstant="${baseDirConstant}/$1"
        fi

        [[ ! -d "$baseDir/faMesh" ]] && mkdir "$baseDir/faMesh" > /dev/null 2>&1
        cp -r "$baseDirConstant/faMesh.org/"* "$baseDir/faMesh" > /dev/null 2>&1

    done
}


####


moveZeroToConstantPolyMesh ()
{
    baseDirZero='0'
    baseDirConstant='constant'
    if [[ ! -z "$1" ]]; then
        baseDirZero="${baseDirZero}/$1"
        baseDirConstant="${baseDirConstant}/$1"
    fi

    [[ ! -d "$baseDirConstant" ]] && mkdir "$baseDirConstant" > /dev/null 2>&1
    mv  "$baseDirZero/polyMesh" "$baseDirConstant/" > /dev/null 2>&1
}

moveZeroToConstantPolyMeshParallel ()
{
    for p in $(ls -1d 'processor'*); do

        baseDirZero="$p/0"
        baseDirConstant="$p/constant"

        if [[ ! -z "$1" ]]; then
            baseDirZero="${baseDirZero}/$1"
            baseDirConstant="${baseDirConstant}/$1"
        fi

        [[ ! -d "$baseDirConstant" ]] && mkdir "$baseDirConstant" > /dev/null 2>&1
        mv  "$baseDirZero/polyMesh" "$baseDirConstant/" > /dev/null 2>&1

    done
}


####


addZeroOrg ()
{
    baseDirZero='0'
    baseDirZeroOrg='0.org'
    if [[ ! -z "$1" ]]; then
        baseDirZero="${baseDirZero}/$1"
        baseDirZeroOrg="${baseDirZeroOrg}/$1"
    fi

    [[ ! -d "$baseDirZero" ]] && mkdir "$baseDirZero" > /dev/null 2>&1

    for f in $(ls -1d "$baseDirZeroOrg/"* 2> /dev/null); do
        if [[ ! -d "$f" ]]; then
            cp "$f" "$baseDirZero" > /dev/null 2>&1
        fi
    done
}

addZeroOrgParallel ()
{
    for p in $(ls -1d 'processor'*); do

        baseDirZero="$p/0"
        baseDirZeroOrg='0.org'

        if [[ ! -z "$1" ]]; then
            baseDirZero="${baseDirZero}/$1"
            baseDirZeroOrg="${baseDirZeroOrg}/$1"
        fi

        [[ ! -d "$baseDirZero" ]] && mkdir "$baseDirZero" > /dev/null 2>&1

        for f in $(ls -1d "$baseDirZeroOrg/"* 2> /dev/null); do
            if [[ ! -d "$f" ]]; then
                cp "$f" "$baseDirZero/." > /dev/null 2>&1
            fi
        done
    done
}


#####


backupToPolyMeshOrg ()
{
    baseDir='constant'
    [[ ! -z "$1" ]] && baseDir="${baseDir}/$1"

    for f in $@; do
        cp "$baseDir/polyMesh/$f" "$baseDir/polyMesh.org" > /dev/null 2>&1
    done
}

backupToFaMeshOrg ()
{
    baseDir='constant'
    [[ ! -z "$1" ]] && baseDir="${baseDir}/$1"

    for f in $@; do
        cp "$baseDir/faMesh/$f" "$baseDir/faMesh.org" > /dev/null 2>&1
    done
}

backupToZeroOrg ()
{
    for f in $@; do
        cp "0/$f" '0.org' > /dev/null 2>&1
    done
}
