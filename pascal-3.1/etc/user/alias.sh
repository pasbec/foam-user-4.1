#!/bin/bash

###############################################################################
PSD="$(dirname "${BASH_SOURCE[0]}")";
sourcePSD () { local psd="$PSD"; source "$1"; PSD="$psd"; }
###############################################################################

alias foamCdUserApps="cd $WM_PROJECT_USER_DIR/applications"
alias foamCdUserRun="cd $WM_PROJECT_USER_DIR/run"
alias foamCdUserTut="cd $WM_PROJECT_USER_DIR/tutorials"
alias foamCdUserSrc="cd $WM_PROJECT_USER_DIR/src"

alias foamCdApps="cd $WM_PROJECT_DIR/applications"
alias foamCdRun="cd $WM_PROJECT_DIR/run"
alias foamCdTut="cd $WM_PROJECT_DIR/tutorials"
alias foamCdSrc="cd $WM_PROJECT_DIR/src"

###############################################################################
