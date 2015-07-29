#!/bin/bash
# July 2015
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Sourcing -------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

PSD="$(dirname "${BASH_SOURCE[0]}")";
sourcePSD () { local psd="$PSD"; source "$1"; PSD="$psd"; }

# --------------------------------------------------------------------------- #

sourcePSD "$PSD/fct_error.sh"

# --------------------------------------------------------------------------- #
# --- Function definitions -------------------------------------------------- #
# --------------------------------------------------------------------------- #

path_exists ()
{
    local path="$1"

    error_fatalIfEmptyVar 'path' "$FUNCNAME" || return 1

    [ -e "$path" ] || return 1

    return 0
}


path_notExists ()
{
    local path="$1"

    error_fatalIfEmptyVar 'path' "$FUNCNAME" || return 1

    [ ! -e "$path" ] || return 1

    return 0
}


path_isAbs ()
{
    local path="$1"

    error_fatalIfEmptyVar 'path' "$FUNCNAME" || return 1

    [ -z "${path##'/'*}" ] || return 1

    return 0
}


path_isNotAbs ()
{
    local path="$1"

    error_fatalIfEmptyVar 'path' "$FUNCNAME" || return 1

    [ ! -z "${path##'/'*}" ] || return 1

    return 0
}


path_normalDir ()
{
    local dir="$1"

    error_fatalIfEmptyVar 'dir' "$FUNCNAME" || return 1

    path_exists "$dir" || return 1

    local normalDir
        normalDir=$(cd "$dir" && pwd) || return 1

    echo "$normalDir"

    return 0
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
