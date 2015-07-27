#!/bin/bash
DIR="$(dirname "${BASH_SOURCE[0]}")"
###############################################################################

source "$DIR/fct_dep.sh"

###############################################################################

alias ll='ls -lAh --color=always --group-directories-first'
alias lL="ls -LlAh --color=always --group-directories-first"

###############################################################################

dep_isHashed 'kate' && alias katen='kate -n'

###############################################################################
