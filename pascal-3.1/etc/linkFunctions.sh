#!/bin/bash

# Absolute paths with python
abspath()
{
  if [[ $(which python) == '' ]]; then
    echo >&2 'ERROR: Python binary not found!'
    return 1
  else
    python -c "import os.path; print os.path.abspath('$1')"
  fi
}

# Relative paths with python
relpath()
{
  if [[ $(which python) == '' ]]; then
    echo >&2 'ERROR: Python binary not found!'
    return 1
  else
    python -c "import os.path; print os.path.relpath('$1','${2:-.}')"
  fi
}

# Translation of a link from absolute to relative form
lnAbsToRel()
{
  local l="$(readlink $1)"
  local dirfile="$(dirname $1)"
  local dirlink="$(dirname $l)"

  if [[ "${l##/*}" == '' ]]; then
    if [[ -e "$l" ]]; then
      local rl="$(relpath $dirlink $dirfile)/$(basename $l)"
      echo "Relinking '$1' -> '${rl#./}'"
      ln -sf "${rl#./}" "$1"
    else
      echo >&2 "ERROR: Absolute link destination $l not found!"
      return 1
    fi
  else
    if [[ -e "$dirfile/$l" ]]; then
      echo "Skipping  '$1'"
    else
      echo >&2 "ERROR: Relative link destination $dirfile/$l not found!"
      return 1
    fi
  fi
}
