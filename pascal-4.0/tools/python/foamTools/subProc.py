#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Module to get easy access to shell subprocesses
# March 2015
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

import os, sys

__path__ = os.path.realpath(__file__)
__dir__ = os.path.dirname(__path__)

import subprocess

# --------------------------------------------------------------------------- #
# --- Function definitions -------------------------------------------------- #
# --------------------------------------------------------------------------- #

def run(exe, arg="", log="", vrb=0):

    # exe : executable string which needs to be hashed
    # arg : arguments string
    # log : log file
    # vrb : verbose level

    # Start sub process and wait till finished
    cmd = exe + " " + arg
    # TODO: Check if "".split() is good enough or not -> shlex.split() ?
    sp = subprocess.Popen(cmd.split())
    #sp = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # FIXME: Fix all stdout/stderr related stuff
    ## Are we verbose?
    #if (vrb > 1):

        ## FIXME: Bug if both stdout and stderr are read simultaneously
        #stderr = ""

        ## Reset lists
        #stdoutL = []
        ##stderrL = []

        ## Poll subprocess and write to console in real-time
        #while sp.poll() is None:

            ## stdout
            #lo = sp.stdout.readline()
            #stdoutL.append(lo)
            #sys.stdout.write(lo)

            ## stderr
            ##le = sp.stderr.readline()
            ##stderrL.append(le)
            ##sys.stderr.write(le)

        ## Save stdout and stderr from list
        #stdout = "".join(stdoutL)
        ##stderr = "".join(stderrL)

    #else:
    if True:

        # Save stdout and stderr via communicate
        stdout, stderr = sp.communicate()

    ## Write log files
    #if not log=="":
        #with open(log + ".out", "w") as logOut:
            #logOut.write(stdout)
        #with open(log + ".err", "w") as logErr:
            #logErr.write(stderr)

    # Check return code and exit if failed
    if not (sp.returncode == 0):
        print(os.path.basename(sys.argv[0]) + ": error: " + "subprocess: \""
              + cmd + "\" returned error code " + str(sp.returncode))
        sys.exit(1)

# --------------------------------------------------------------------------- #
# --- End of module --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

