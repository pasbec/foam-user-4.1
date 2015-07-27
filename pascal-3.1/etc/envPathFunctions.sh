#!/bin/bash

# Environment variable tools
env_set ()
{
    local ENV="$1"
    local KEY="$2"
    eval "export $ENV='$KEY'"
}

env_unset ()
{
    local ENV="$1"
    unset "$ENV"
}

env_unset_regex ()
{
    local KEY_REGEX="$1"
    local KEY_EXPAND="$(env | sed -n "/$KEY_REGEX=.*/p" | sed 's/^\([^=]*\)=\(.*\)/\1/g')"
    for var in $KEY_EXPAND; do
        unset "$var"
    done
}

env_append ()
{
    local ENV="$1"
    local KEY="$2"
    env_remove "$ENV" "$KEY"
    if [[ -z "$(eval "echo \$$ENV")" ]]; then
        eval "export $ENV=\"$KEY\""
    else
        eval "export $ENV=\"\$$ENV:$KEY\""
    fi
}
env_prepend ()
{
    local ENV="$1"
    local KEY="$2"
    env_remove "$ENV" "$KEY"
    if [[ -z "$(eval "echo \$$ENV")" ]]; then
        eval "export $ENV=\"$KEY\""
    else
        eval "export $ENV=\"$KEY:\$$ENV\""
    fi
}
env_remove ()
{
    local ENV="$1"
    local KEY="$2"
    eval "$ENV=\$(echo -n \$$ENV | awk -v RS=: -v ORS=: '\$0 != \"$KEY\"' | sed 's/:$//')"
    [[ -z "$(eval "echo \$$ENV")" ]] && unset "$ENV"
}
env_removeany ()
{
    local ENV="$1"
    local KEY="$2"
    eval "$ENV=\$(echo -n \$$ENV | sed 's/^[^:]*$KEY[^:]*://g' | sed 's/[^:]*$KEY[^:]*://g' | sed 's/[^:]*$KEY[^:]*$//g')"
    [[ -z "$(eval "echo \$$ENV")" ]] && unset "$ENV"
}
env_removeall ()
{
    local KEY="$1"
    for var in $(env | grep "$KEY"); do
        local ENV="$(echo $var | sed 's/^\([^=]*\)=\(.*\)/\1/g')"
        env_removeany "$ENV" "$KEY"
    done
    for var in $(env | grep '=$'); do
        local ENV="$(echo $var | sed 's/^\([^=]*\)=\(.*\)/\1/g')"
        env_unset "$ENV"
    done
}

# Path variable tools
path_append ()  { env_append 'PATH' "$1"; }
path_prepend ()  { env_prepend 'PATH' "$1"; }
path_remove ()  { path_remove 'PATH' "$1"; }
path_removeany ()  { path_removeany 'PATH' "$1"; }

# Library path variable tools
lpath_append ()  { env_append 'LIBRARY_PATH' "$1"; }
lpath_prepend ()  { env_prepend 'LIBRARY_PATH' "$1"; }
lpath_remove ()  { path_remove 'LIBRARY_PATH' "$1"; }
lpath_removeany ()  { path_removeany 'LIBRARY_PATH' "$1"; }

# LD Library path variable tools
ldpath_append ()  { env_append 'LD_LIBRARY_PATH' "$1"; }
ldpath_prepend ()  { env_prepend 'LD_LIBRARY_PATH' "$1"; }
ldpath_remove ()  { path_remove 'LD_LIBRARY_PATH' "$1"; }
ldpath_removeany ()  { path_removeany 'LD_LIBRARY_PATH' "$1"; }

# Python path variable tools
pypath_append ()  { env_append 'PYTHONPATH' "$1"; }
pypath_prepend ()  { env_prepend 'PYTHONPATH' "$1"; }
pypath_remove ()  { pypath_remove 'PYTHONPATH' "$1"; }
pypath_removeany ()  { pypath_removeany 'PYTHONPATH' "$1"; }

# C include path
cipath_append ()  { env_append 'C_INCLUDE_PATH' "$1"; }
cipath_prepend ()  { env_prepend 'C_INCLUDE_PATH' "$1"; }
cipath_remove ()  { pipath_remove 'C_INCLUDE_PATH' "$1"; }
cipath_removeany ()  { pipath_removeany 'C_INCLUDE_PATH' "$1"; }

# C++ include path
cppipath_append ()  { env_append 'CPLUS_INCLUDE_PATH' "$1"; }
cppipath_prepend ()  { env_prepend 'CPLUS_INCLUDE_PATH' "$1"; }
cppipath_remove ()  { pipath_remove 'CPLUS_INCLUDE_PATH' "$1"; }
cppipath_removeany ()  { pipath_removeany 'CPLUS_INCLUDE_PATH' "$1"; }