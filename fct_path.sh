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

    error_fatalIfEmptyVar 'path' "$FUNCNAME" || return $?

    [ -e "$path" ] || return $?

    return 0
}


path_notExists ()
{
    local path="$1"

    error_fatalIfEmptyVar 'path' "$FUNCNAME" || return $?

    [ ! -e "$path" ] || return $?

    return 0
}


path_isAbs ()
{
    local path="$1"

    error_fatalIfEmptyVar 'path' "$FUNCNAME" || return $?

    [ -z "${path##'/'*}" ] || return $?

    return 0
}


path_isNotAbs ()
{
    local path="$1"

    error_fatalIfEmptyVar 'path' "$FUNCNAME" || return $?

    [ ! -z "${path##'/'*}" ] || return $?

    return 0
}


path_normalDir ()
{
    local dir="$1"

    error_fatalIfEmptyVar 'dir' "$FUNCNAME" || return $?

    path_exists "$dir" || return $?

    local normalDir
        normalDir=$(cd "$dir" && pwd) || return $?

    echo "$normalDir" || return $?

    return 0
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
