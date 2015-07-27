#!/bin/bash
DIR="$(dirname "${BASH_SOURCE[0]}")"
###############################################################################

source "$DIR/fct_error.sh"

###############################################################################

bool_unset ()
{
    local name="$1"
    local prefix="${2:-bool}"

    error_fatalIfEmptyVar 'name prefix' "$FUNCNAME" || return 1

    eval "unset -f ${prefix}_${name}" || return 1

    return 0
}


bool_set ()
{
    local name="$1"
    local value="$2"
    local prefix="${3:-bool}"

    error_fatalIfEmptyVar 'name value prefix' "$FUNCNAME" || return 1

    if [ $value -eq 0 ]; then
        eval "${prefix}_${name} () { return 1; }" || return 1

    elif [ $value -eq 1 ]; then
        eval "${prefix}_${name} () { return 0; }" || return 1
    else
        error_fatal \
            "Boolean value may only take '0' or '1' (given: '$value')." \
            "$FUNCNAME" || return 1
    fi

    return 0
}


bool_true ()
{
    local name="$1"
    local prefix="${2:-bool}"

    error_fatalIfEmptyVar 'name prefix' "$FUNCNAME" || return 1

    bool_set "$name" 1 "$prefix" || return 1

    return 0
}


bool_false ()
{
    local name="$1"
    local prefix="${2:-bool}"

    error_fatalIfEmptyVar 'name prefix' "$FUNCNAME" || return 1

    bool_set "$name" 0 "$prefix" || return 1

    return 0
}

###############################################################################
