#!/usr/bin/env python

import os
import sys
import pylauncher

cmdLineFile = sys.argv[1]
numCoresPerCmd = 4
if len(sys.argv) == 3:
    numCoresPerCmd = int(sys.argv[2])

pylauncher.ClassicLauncher( cmdLineFile, cores=numCoresPerCmd )

sys.exit()
