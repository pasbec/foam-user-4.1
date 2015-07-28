#!/bin/bash

# TODO [High]: Check if flag exists!

###############################################################################
DIR="$(dirname "${BASH_SOURCE[0]}")"
###############################################################################

source "$DIR/fct_error.sh"
source "$DIR/fct_bool.sh"

###############################################################################
DIR="$(dirname "${BASH_SOURCE[0]}")"
###############################################################################

flag_nameOkOrExit ()
{
    local name="$1"

    error_fatalIfEmptyVar 'name' "$FUNCNAME" || return 1

    local nameStripped="$name"

    nameStripped=${name%%*'_'}
    nameStripped=${name%%'_'*}
    nameStripped=${name##*'_'}
    nameStripped=${name##'_'*}

    [ "$name" = "$nameStripped" ] || \
        error_fatal "Flag name $name contains 'underscore(s)'." "$FUNCNAME"

    [ $? -eq 0 ] || return 1

    return 0
}


flag_init ()
{
    local flags="$1"

    error_fatalIfEmptyVar 'flags' "$FUNCNAME" || return 1

    eval "flag_flagList='$flags'" || return 1

    local f
    for f in $flags; do
        flag_nameOkOrExit "$f" || return 1
        bool_false "$f" 'flag_isSet' || return 1
    done

    return 0
}

flag_initRead ()
{
    local flags="$1"
    shift

    error_fatalIfEmptyVar 'flags' "$FUNCNAME" || return 1

    flag_init "$flags"

    while [ $# -gt 0 ]; do
        flag_enable "$1"
        shift
    done

    return 0
}


flag_set ()
{
    local name="$1"
    local value="$2"

    error_fatalIfEmptyVar 'name value' "$FUNCNAME" || return 1

    local f
    for f in $flag_flagList; do

        if [ "$f" = "$name" ]; then

	    bool_set "$f" "$value" 'flag_isSet' || return 1
            return 0
        fi
    done

    error_fatal "Flag '$name' is unknown." "$FUNCNAME" || return 1
}


flag_enable ()
{
    local name="$1"

    error_fatalIfEmptyVar 'name' "$FUNCNAME" || return 1

    flag_set "$name" 1 || return 1

    return 0
}


flag_disable ()
{
    local name="$1"

    error_fatalIfEmptyVar 'name' "$FUNCNAME" || return 1

    flag_set "$name" 0 || return 1

    return 0
}

###############################################################################
