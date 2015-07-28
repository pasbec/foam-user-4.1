#!/bin/bash

###############################################################################
PSD="$(dirname "${BASH_SOURCE[0]}")";
sourcePSD () { local psd="$PSD"; source "$1"; PSD="$psd"; }
###############################################################################

for fct in $(ls -1 "$PSD/fct_"* 2>/dev/null); do
    sourcePSD "$fct"
done

sourcePSD "$PSD/colors.sh"
sourcePSD "$PSD/alias.sh"

###############################################################################
