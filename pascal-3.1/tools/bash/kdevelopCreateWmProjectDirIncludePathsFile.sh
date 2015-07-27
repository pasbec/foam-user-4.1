#!/bin/bash

# Created/Updated file .kdev_include_paths in directory $WM_PROJECT_DIR/src for KDevelop

if [ -w $WM_PROJECT_DIR ]; then
  find "$WM_PROJECT_DIR/src" -name 'lnInclude' -print \
      | tee "$WM_PROJECT_DIR/src/.kdev_include_paths" \
      | tee "$WM_PROJECT_DIR/applications/.kdev_include_paths" \
      > /dev/null
  chmod go+r "$WM_PROJECT_DIR/src/.kdev_include_paths"
  chmod go+r "$WM_PROJECT_DIR/applications/.kdev_include_paths"
  echo "Created/Updated file .kdev_include_paths in directory $WM_PROJECT_DIR/src and $WM_PROJECT_DIR/applications."
  exit 0
else
  echo "You do not have write permissions for \$WM_PROJECT_DIR ($WM_PROJECT_DIR)."
  echo "Try to rerun this script with 'sudo -E'."
  echo "WARNING: Pay attention to the -E (preserve-environment) option!"
  exit 1
fi