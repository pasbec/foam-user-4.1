#!/bin/bash
DIR="$(dirname "${BASH_SOURCE[0]}")"
###############################################################################

env_set COLOR_ESC_SEQ  "\x1b["
env_set COLOR_RESET    $COLOR_ESC_SEQ"39;49;00m"
env_set COLOR_RED      $COLOR_ESC_SEQ"31;01m"
env_set COLOR_GREEN    $COLOR_ESC_SEQ"32;01m"
env_set COLOR_YELLOW   $COLOR_ESC_SEQ"33;01m"
env_set COLOR_BLUE     $COLOR_ESC_SEQ"34;01m"
env_set COLOR_MAGENTA  $COLOR_ESC_SEQ"35;01m"
env_set COLOR_CYAN     $COLOR_ESC_SEQ"36;01m"

###############################################################################
