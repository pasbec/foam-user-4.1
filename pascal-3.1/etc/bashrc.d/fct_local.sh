#!/bin/bash

###############################################################################
DIR="$(dirname "${BASH_SOURCE[0]}")"
###############################################################################

local_set ()
{
    eval "$1='$2'" || return 1

    return 0
}

###############################################################################
