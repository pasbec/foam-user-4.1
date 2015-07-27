#!/bin/bash

reqCheckList='cpp awk grep mpirun python mayavi2 sed'

# Check gloabl etc directory
[[ ! -d "$globalEtcDir" ]] && { echo >&2 "ERROR: The global etc directory $globalEtcDir is not existing."; return 1; }

# Check for necessary apps
checkBin () {
  [[ -z "$1" ]] && { echo >&2 "ERROR: No binary name given to check."; return 1; }
  hash $1 2>/dev/null || { echo >&2 "WARNING: The application/tool $1 is not available."; return 1; }
}

for req in $reqCheckList; do
  checkBin $req
done
