#!/usr/bin/python
# Init script for foamToolLib
# March 2015
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

# pyFoam module paths
# FIXME: Make kdevelop syntax highlighting possible without $PYTHONPATH env
# sys.path.append(os.environ['WM_PROJECT_USER_DIR'] + '/tools/python/pyFoam/lib/python2.7/site-packages')

from foamTools import fileEdit
from foamTools import subProc
from foamTools import ioInfo
from foamTools import classInfo
from foamTools import dictHelper
from foamTools import fieldWriteFromRaw

# --------------------------------------------------------------------------- #
# --- End of init ----------------------------------------------------------- #
# --------------------------------------------------------------------------- #

