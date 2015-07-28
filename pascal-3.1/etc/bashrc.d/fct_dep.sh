#!/bin/bash

###############################################################################
PSD="$(dirname "${BASH_SOURCE[0]}")";
sourcePSD () { local psd="$PSD"; source "$1"; PSD="$psd"; }
###############################################################################

sourcePSD "$PSD/fct_error.sh"

###############################################################################

dep_isHashed ()
{
    local bin="$1"

    error_fatalIfEmptyVar 'bin' "$FUNCNAME" || return 1

    hash $bin 2>/dev/null || return 1

    return 0
}


dep_isHashedVerb ()
{
    local bin="$1"

    error_fatalIfEmptyVar 'bin' "$FUNCNAME" || return 1

    dep_isHashed $bin || \
        error_warning "The binary/application '$bin' is not hashed." \
        "$FUNCNAME" || return 1

    return $?
}


dep_isHashedOrExit ()
{
    local bin="$1"

    error_fatalIfEmptyVar 'bin' "$FUNCNAME" || return 1

    local msg
        msg="This is considered fatal."
        msg="$msg Please make sure to install '$bin'."

    dep_isHashedVerb "$bin" || error_fatal "$msg" "$FUNCNAME" || return 1

    return 0
}


dep_listIsHashed ()
{
    local binList="$1"

    error_fatalIfEmptyVar 'binList' "$FUNCNAME" || return 1

    local err=0

    local b
    for b in $binList; do
        dep_isHashed $b
        [ $? -eq 0 ] || err=1
    done

    return $err
}


dep_listIsHashedVerb ()
{
    local binList="$1"

    error_fatalIfEmptyVar 'binList' "$FUNCNAME" || return 1

    local err=0

    local b
    for b in $binList; do
        dep_isHashedVerb $b
        [ $? -eq 0 ] || err=1
    done

    return $err
}


dep_listIsHashedOrExit ()
{
    local binList="$1"

    error_fatalIfEmptyVar 'binList' "$FUNCNAME" || return 1

    local b
    for b in $binList; do
        dep_isHashedOrExit $b || return 1
    done

    return 0
}

###############################################################################
