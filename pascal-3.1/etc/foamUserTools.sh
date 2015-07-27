#!/bin/bash

foamUserToolBinLinksRemoveDangling() {
  for l in $(find "$WM_PROJECT_USER_DIR/bin" -type l); do
    if [[ ! -e "$(readlink $(basename $l))" ]]; then
      rm $l
    fi
  done
}

# Script links for binary creation
foamUserToolBinLinks() {
  local scriptdir="$WM_PROJECT_USER_DIR/tools"
  local bindir="$WM_PROJECT_USER_DIR/bin"
  local progtag='foamUserTool'

  foamUserToolBinLinksRemoveDangling

  for file in $(find $scriptdir -maxdepth 2 -type f 2> /dev/null); do
    if [[ -x $file ]]; then
      local fileBase="$(basename $file)"
      local source="$file"
      local fileNoExt="${fileBase%.*}"
      local target="$bindir/${progtag}${fileNoExt^}"
      echo "Linking $file to $target"
      if [[ -h "$target" ]]; then
        rm "$target"
      elif [[ -e "$target" ]]; then
        echo >&2 "ERROR: Link target $target exists and is no link!"
        return 1
      fi
      ln -s "$source" "$target"
      lnAbsToRel "$target"
      if [[ $? != 0 ]]; then
        echo >&2 'ERROR: Linking failed!'
        return 1
      fi
    else
      echo "Skipping $file"
    fi
  done
  return 0
}

# Link tool scripts to binary folder and add the latter to path
foamUserToolBinLinks 1> /dev/null || return $?
env_prepend 'PATH' "$WM_PROJECT_USER_DIR/bin"

# Add python library (module) path to PYTHONPATH variable
env_prepend 'PYTHONPATH' "$WM_PROJECT_USER_DIR/tools/python"
