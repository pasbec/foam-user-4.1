#!/bin/bash

###############################################################################
src () { [ ! -z "$1" ] && . "$1"; PSD="$(dirname "${BASH_SOURCE[0]}")"; }; src
###############################################################################

src "$PSD/fct_error.sh"
src "$PSD/fct_dep.sh"

###############################################################################

env_set ()
{
    local var="$1"
    local key="$2"

    error_fatalIfEmptyVar 'var key' "$FUNCNAME" || return 1

    eval "export $var='$key'" || return 1

    return 0
}


env_unset ()
{
    local var="$1"

    error_fatalIfEmptyVar 'var' "$FUNCNAME" || return 1

    unset "$var" || return 1

    return 0
}


env_unsetRegex ()
{
    local key_regex="$1"

    error_fatalIfEmptyVar 'key_regex' "$FUNCNAME" || return 1

    if dep_isHashedOrExit 'sed'; then

        local key_expand
            key_expand="$( env | \
            sed -n "/$key_regex=.*/p" | \
            sed 's/^\([^=]*\)=\(.*\)/\1/g' )" || return 1

        local v
        for v in $key_expand; do
            unset "$v" || return 1
        done

    else
        return 1
    fi

    return 0
}


env_remove ()
{
    local var="$1"
    local key="$2"

    error_fatalIfEmptyVar 'var key' "$FUNCNAME" || return 1

    if dep_listIsHashedOrExit 'sed awk'; then

# TODO [High]: Error handling and formatting
        eval "$var=\$(echo -n \$$var | awk -v RS=: -v ORS=: '\$0 != \"$key\"' | sed 's/:$//')"
        [ -z "$(eval "echo \$$var")" ] && unset "$var"

    else
        return 1
    fi

    return 0
}


env_removeAny ()
{
    local var="$1"
    local key="$2"

    error_fatalIfEmptyVar 'var key' "$FUNCNAME" || return 1

    if dep_isHashedOrExit 'sed'; then

# TODO [High]: Error handling and formatting
        eval "$var=\$(echo -n \$$var | sed 's/^[^:]*$key[^:]*://g' | sed 's/[^:]*$key[^:]*://g' | sed 's/[^:]*$key[^:]*$//g')"
        [ -z "$(eval "echo \$$var")" ] && unset "$var"

    else
        return 1
    fi

    return 0
}


env_removeAll ()
{
    local key="$1"

    error_fatalIfEmptyVar 'key' "$FUNCNAME" || return 1

    if dep_listIsHashedOrExit 'grep sed'; then

        local varList
            varList="$(env | grep "$key" | sed 's/^\([^=]*\)=\(.*\)/\1/g')" \
            || return 1

        local v
        for v in $varList; do

            env_removeAny "$v" "$key" || return 1

        done

        local varList
            varList="$(env | grep '=$' | sed 's/=//g')" \
            || return 1

        local v
        for v in $varList; do

            env_unset "$v" || return 1

        done

    else
        return 1
    fi

    return 0
}


env_append ()
{
    local var="$1"
    local key="$2"

    error_fatalIfEmptyVar 'var key' "$FUNCNAME" || return 1

    env_remove "$var" "$key" || return 1

# TODO [High]: Error handling and formatting
    if [ -z "$(eval "echo \$$var")" ]; then
        eval "export $var=\"$key\""
    else
        eval "export $var=\"\$$var:$key\""
    fi

    return 0
}


env_prepend ()
{
    local var="$1"
    local key="$2"

    error_fatalIfEmptyVar 'var key' "$FUNCNAME" || return 1

    env_remove "$var" "$key" || return 1

# TODO [High]: Error handling and formatting
    if [ -z "$(eval "echo \$$var")" ]; then
        eval "export $var=\"$key\""
    else
        eval "export $var=\"$key:\$$var\""
    fi

    return 0
}
