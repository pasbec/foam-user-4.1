#!/bin/bash

###############################################################################
PSD="$(dirname "${BASH_SOURCE[0]}")";
sourcePSD () { local psd="$PSD"; source "$1"; PSD="$psd"; }
###############################################################################

sourcePSD "$PSD/env.sh"
sourcePSD "$PSD/dep.sh"
sourcePSD "$PSD/alias.sh"

###############################################################################