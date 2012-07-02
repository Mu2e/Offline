#!/bin/bash
#
# Check the status of packages of interest
#
# $Id: update.sh,v 1.3 2012/07/02 15:50:50 kutschke Exp $
# $Date: 2012/07/02 15:50:50 $
# $Author: kutschke $
#
# Contact person Rob Kutschke
#

for PKG in BField BaBar BbrGeom CLHEP DetectorModel KalmanTrack MatEnv ProbTools TrajGeom TrkBase difAlgebra
do
  echo "Checking for updates in " $PKG
  svn update BaBar/$PKG
done
