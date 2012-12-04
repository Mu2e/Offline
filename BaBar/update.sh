#!/bin/bash
#
# Check the status of packages of interest
#
# $Id: update.sh,v 1.5 2012/12/04 00:51:28 tassiell Exp $
# $Date: 2012/12/04 00:51:28 $
# $Author: tassiell $
#
# Contact person Rob Kutschke
#

if ! [[ -f "SConstruct" ]]; then
  echo "ERROR: must be executed from top level of offline code." 1>&2
  exit 1
fi

for PKG in BField BaBar BbrGeom CLHEP DetectorModel KalmanTrack MatEnv ProbTools TrajGeom TrkBase difAlgebra Dch
do
  echo "Checking for updates in $PKG"
  svn update "BaBar/$PKG"
done
