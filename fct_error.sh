#!/bin/bash
DIR="$(dirname "${BASH_SOURCE[0]}")"
###############################################################################

error_fatal ()
{
    local message="$1"
    local function="$2"

    local functionStr=''
    [ ! -z "$function" ] && functionStr="[$function] "

    echo >&2 "ERROR: $functionStr$message"

    case $- in
        *i*)
            return 1
        ;;
        *)
            exit 1
        ;;
    esac
}


error_fatalIfEmptyVar ()
{
    local variable="$1"
    local function="$2"

    local v
    for v in $variable; do

        local value
            value="$(eval "echo \$$v")" || return 1

        [ ! -z "$value" ] || \
            error_fatal "Variable '$v' is empty." "$function"

        [ $? -eq 0 ] || return 1

    done

    return 0
}


error_warning ()
{
    local message="$1"
    local function="$2"

    local functionStr=''
    [ ! -z "$function" ] && functionStr="[$function] "

    echo >&2 "WARNING: $functionStr$message"

    return 1
}


error_info ()
{
    local message="$1"
    local function="$2"

    local functionStr=''
    [ ! -z "$function" ] && functionStr="[$function] "

    echo >&2 "INFO: $functionStr$message"

    return 0
}


error_echo ()
{
    local message="$1"
    local function="$2"

    local functionStr=''
    [[ ! -z "$function" ]] && functionStr="[$function] "

    echo >&2 "$functionStr$message"

    return 0
}

###############################################################################
