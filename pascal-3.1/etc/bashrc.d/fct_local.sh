#!/bin/bash

###############################################################################
DIR="$(dirname "${BASH_SOURCE[0]}")"
###############################################################################

local_sourceDir ()
{
    local var="$1"

    eval "$var="$(dirname "${BASH_SOURCE[0]}")""

    return $?
}

local_source ()
{
    local file="$1"
    local var="$2"

    source "$file" || return 1

    local_sourceDir "$var"

    return $?
}

local_set ()
{
    eval "$1='$2'"

    return $?
}

###############################################################################
