#!/bin/bash

# Basic set of color variables
env_set COL_ESC_SEQ  "\x1b["
env_set COL_RESET    $COL_ESC_SEQ"39;49;00m"
env_set COL_RED      $COL_ESC_SEQ"31;01m"
env_set COL_GREEN    $COL_ESC_SEQ"32;01m"
env_set COL_YELLOW   $COL_ESC_SEQ"33;01m"
env_set COL_BLUE     $COL_ESC_SEQ"34;01m"
env_set COL_MAGENTA  $COL_ESC_SEQ"35;01m"
env_set COL_CYAN     $COL_ESC_SEQ"36;01m"
