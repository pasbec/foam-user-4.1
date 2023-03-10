#!/bin/bash
# July 2015
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Sourcing -------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

PSD="$(dirname "${BASH_SOURCE[0]}")";
sourcePSD () { local psd="$PSD"; source "$1"; PSD="$psd"; }

# --------------------------------------------------------------------------- #

sourcePSD "$PSD/fct_local.sh"

# --------------------------------------------------------------------------- #
# --- Script ---------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

local_set COLOR_ESC_SEQ  "\x1b["
local_set COLOR_RESET    $COLOR_ESC_SEQ"39;49;00m"
local_set COLOR_RED      $COLOR_ESC_SEQ"31;01m"
local_set COLOR_GREEN    $COLOR_ESC_SEQ"32;01m"
local_set COLOR_YELLOW   $COLOR_ESC_SEQ"33;01m"
local_set COLOR_BLUE     $COLOR_ESC_SEQ"34;01m"
local_set COLOR_MAGENTA  $COLOR_ESC_SEQ"35;01m"
local_set COLOR_CYAN     $COLOR_ESC_SEQ"36;01m"

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
